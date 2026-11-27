/*
  CRAM move prior for Stockfish (experimental add-on, not part of upstream Stockfish)

  A CRAM net is a small, sparse "classical energy" over (position, move): for every legal
  move the score is a sum of pairwise terms between the move's from/to squares and the
  pieces on the board, plus a few attack-count terms. softmax(score) over the legal moves
  is a cheap move prior ("policy"). Training is done offline with chess_cram.cpp
  (pseudolikelihood on Stockfish MultiPV labels); this file only loads and evaluates a net.

  The net is 8-bit quantized: int8 weights plus one integer fixed-point multiplier per parameter
  group. Inference uses INTEGER arithmetic only (int8 weights summed in int32, scaled by the
  integer multiplier into int64 fixed-point scores, integer exp table for the gate), so results
  are bit-exact on every platform and compiler: no floating point, no rounding effects.
  The binary layout below MUST stay identical to chess_cram.cpp.

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/

#ifndef CRAM_H_INCLUDED
#define CRAM_H_INCLUDED

#include <optional>
#include <string>

#include "movegen.h"
#include "types.h"

namespace Stockfish {

class Position;

namespace CRAM {

// File format "CRM3" (little-endian): u32 magic, u32 NPARAM, u32 NGROUPS, u32 FRAC,
// i32 mult[NGROUPS], int8 w[NPARAM]. The weight w[i] of group g has the value w[i] * mult[g] / 2^FRAC.
// A move score is sum_g mult[g] * (integer sum of the int8 weights of group g), in units of 2^-FRAC.
constexpr u32 MAGIC     = 0x43524d33;  // "3MRC": bytes 33 4d 52 43
constexpr u32 MAGIC_Q8F = 0x43524d32;  // old int8 + float-scale format, recognized only to give an error
constexpr u32 MAGIC_F32 = 0x43524d31;  // old float format, recognized only to give an error
constexpr int FRAC      = 20;          // fixed-point fraction bits of scores and group multipliers
constexpr i32 MAX_MULT  = 1 << 24;

// Name of the default net. It is embedded into the binary at build time (see cram.cpp).
constexpr const char* DefaultFile = "model-i8-quiets2.bin";

constexpr int NF      = 768;  // 12 piece kinds (side-to-move relative) x 64 squares
constexpr int OFF_WF  = 0;
constexpr int OFF_WT  = OFF_WF + 64 * NF;
constexpr int OFF_K   = OFF_WT + 64 * NF;
constexpr int OFF_BF  = OFF_K + 4096;
constexpr int OFF_BT  = OFF_BF + 64;
constexpr int OFF_P   = OFF_BT + 64;
constexpr int OFF_KMC = OFF_P + 5;
constexpr int OFF_A   = OFF_KMC + 6 * 7;
constexpr int OFF_AF  = OFF_A + 6 * 4 * 4;
constexpr int NPARAM  = OFF_AF + 6 * 4;

// Parameter groups, each with its own scale: [GOFF[g], GOFF[g+1]) are the weights of group g.
constexpr int NGROUPS = 9;
constexpr int GOFF[NGROUPS + 1] = {OFF_WF, OFF_WT, OFF_K, OFF_BF, OFF_BT, OFF_P, OFF_KMC, OFF_A, OFF_AF, NPARAM};
enum Group { G_WF, G_WT, G_K, G_BF, G_BT, G_P, G_KMC, G_A, G_AF };

// Loading / configuration. These must only be called while no search is running
// (the UCI option callbacks wait for the search to finish first).
// Loads a net from a UTF-8 path (or the embedded net for DefaultFile); returns an info message
std::optional<std::string> load(const std::string& path);
void                       unload();
bool                       loaded();

void set_weight(int w);  // move ordering strength, in history-value units per nat
void set_plies(int p);   // use the prior at plies < p
void set_gate(int pct);  // 0 = off; otherwise play the CRAM move at the root if p_max >= pct %
int  gate_percent();

// Move ordering: true if the prior should be used at this ply.
inline bool ordering_active(int ply, int depth);

// Adds weight * (score - mean score) to the value of every move in [first, last).
void add_ordering_bonus(const Position& pos, ExtMove* first, ExtMove* last, int depth);

// Root gate: scores the candidate moves, returns false if no net is loaded or no candidates.
// On success 'best' is the most likely move and 'confidencePpm' its softmax probability in
// millionths (1000000 = 100 %), computed with integers only. Ties go to the smaller Move::raw().
bool gate_pick(const Position& pos, const Move* cands, int n, Move& best, i64& confidencePpm);

// Text report for the "cram" debug command (all legal moves with their scores).
std::string report(const Position& pos, bool chess960);

// ---- inline part -------------------------------------------------------------------------
extern bool g_loaded;
extern int  g_plies;
inline bool ordering_active(int /*ply*/, int depth) { return g_loaded && depth <= 3; } // && ply < g_plies && depth > ply; }

}  // namespace CRAM
}  // namespace Stockfish

#endif  // #ifndef CRAM_H_INCLUDED
