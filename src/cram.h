/*
  CRAM move prior for Stockfish (experimental add-on, not part of upstream Stockfish)

  A CRAM net is a small, sparse "classical energy" over (position, move): for every legal
  move the score is a sum of pairwise terms between the move's from/to squares and the
  pieces on the board, plus a few attack-count terms. softmax(score) over the legal moves
  is a cheap move prior ("policy"). Training is done offline with chess_cram.cpp
  (pseudolikelihood on Stockfish MultiPV labels); this file only loads and evaluates a net.

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

// File format: u32 magic, u32 parameter count, then NPARAM little-endian IEEE-754 floats.
constexpr u32 MAGIC = 0x43524d31;  // "1MRC" little-endian bytes: 31 4d 52 43

constexpr int NF     = 768;  // 12 piece kinds (side-to-move relative) x 64 squares
constexpr int OFF_WF = 0;
constexpr int OFF_WT = OFF_WF + 64 * NF;
constexpr int OFF_K  = OFF_WT + 64 * NF;
constexpr int OFF_BF = OFF_K + 4096;
constexpr int OFF_BT = OFF_BF + 64;
constexpr int OFF_P  = OFF_BT + 64;
constexpr int OFF_KMC = OFF_P + 5;
constexpr int OFF_A   = OFF_KMC + 6 * 7;
constexpr int OFF_AF  = OFF_A + 6 * 4 * 4;
constexpr int NPARAM  = OFF_AF + 6 * 4;

// Loading / configuration. These must only be called while no search is running
// (the UCI option callbacks wait for the search to finish first).
std::optional<std::string> load(const std::string& path);  // returns an info message
void                       unload();
bool                       loaded();

void set_weight(int w);   // move ordering strength, in history-value units per nat
void set_plies(int p);    // use the prior at plies < p
void set_gate(int pct);   // 0 = off; otherwise play the CRAM move at the root if p_max >= pct %
int  gate_percent();

// Move ordering: true if the prior should be used at this ply.
inline bool ordering_active(int ply);

// Adds weight * (score - mean score) to the value of every move in [first, last).
void add_ordering_bonus(const Position& pos, ExtMove* first, ExtMove* last);

// Root gate: scores the candidate moves, returns false if no net is loaded or no candidates.
// On success 'best' is the most likely move and 'confidence' its softmax probability.
bool gate_pick(const Position& pos, const Move* cands, int n, Move& best, float& confidence);

// Text report for the "cram" debug command (all legal moves with their scores).
std::string report(const Position& pos, bool chess960);

// ---- inline part -------------------------------------------------------------------------
extern bool g_loaded;
extern int  g_plies;
inline bool ordering_active(int ply) { return g_loaded && ply < g_plies; }

}  // namespace CRAM
}  // namespace Stockfish

#endif  // #ifndef CRAM_H_INCLUDED
