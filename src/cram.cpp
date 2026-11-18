/*
  CRAM move prior for Stockfish (experimental add-on, not part of upstream Stockfish)

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/

#include "cram.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

#include "bitboard.h"
#include "movegen.h"
#include "position.h"
#include "uci.h"

namespace Stockfish::CRAM {

bool g_loaded = false;
int  g_plies  = 4;

namespace {

std::vector<float> g_w;  // the net, NPARAM floats; only replaced while no search runs
int                g_weight = 256;
int                g_gate   = 0;

// Per-position evaluation context. Field scores are computed lazily and shared by all moves
// that start/end on the same square (about 2x fewer dot products than moves).
struct Ctx {
    Ctx(const Position& p) :
        pos(p),
        us(p.side_to_move()),
        flip(us == BLACK ? 56 : 0) {
        Bitboard b = pos.pieces();
        while (b)
        {
            Square s  = pop_lsb(b);
            Piece  pc = pos.piece_on(s);
            fidx[nfeat++] =
              u16((type_of(pc) - 1 + 6 * (color_of(pc) != us)) * 64 + (int(s) ^ flip));
        }
    }

    float field(int base, int bias, int o) const {
        const float* w   = g_w.data();
        const float* row = w + base + o * NF;
        float        v   = w[bias + o];
        for (int j = 0; j < nfeat; ++j)
            v += row[fidx[j]];
        return v;
    }

    float score(Move m) {
        const float* w    = g_w.data();
        Square       from = m.from_sq(), to = m.to_sq();
        const bool   castling = m.type_of() == CASTLING;
        if (castling)  // Stockfish encodes castling as "king takes own rook"; use the king's target
            to = make_square(to > from ? FILE_G : FILE_C, rank_of(from));

        const int of = int(from) ^ flip, ot = int(to) ^ flip;
        if (!(doneF >> of & 1))
            fF[of] = field(OFF_WF, OFF_BF, of), doneF |= u64(1) << of;
        if (!(doneT >> ot & 1))
            fT[ot] = field(OFF_WT, OFF_BT, ot), doneT |= u64(1) << ot;

        const int type = type_of(pos.piece_on(from)) - 1;
        int       cap  = 6;
        if (!castling)
        {
            if (m.type_of() == EN_PASSANT)
                cap = 0;
            else if (pos.piece_on(to) != NO_PIECE)
                cap = type_of(pos.piece_on(to)) - 1;
        }
        const int promo = m.type_of() == PROMOTION ? m.promotion_type() - 1 : 0;

        const Bitboard att = pos.attackers_to(to);
        const int      eN  = std::min(3, popcount(att & pos.pieces(~us)));
        const int      oN  = std::min(3, popcount(att & pos.pieces(us)));
        const int      eF  = std::min(3, popcount(pos.attackers_to(from) & pos.pieces(~us)));

        return fF[of] + fT[ot] + w[OFF_K + of * 64 + ot] + (promo ? w[OFF_P + promo] : 0.f)
             + w[OFF_KMC + type * 7 + cap] + w[OFF_A + (type * 4 + eN) * 4 + oN]
             + w[OFF_AF + type * 4 + eF];
    }

    const Position& pos;
    Color           us;
    int             flip;
    int             nfeat = 0;
    u16             fidx[32];
    u64             doneF = 0, doneT = 0;
    float           fF[64], fT[64];
};

}  // namespace

bool loaded() { return g_loaded; }

void unload() {
    g_loaded = false;
    std::vector<float>().swap(g_w);
}

void set_weight(int w) { g_weight = w; }
void set_plies(int p) { g_plies = p; }
void set_gate(int pct) { g_gate = pct; }
int  gate_percent() { return g_gate; }

std::optional<std::string> load(const std::string& path) {
    if (path.empty() || path == "<empty>")
    {
        unload();
        return "CRAM net unloaded";
    }

    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f)
        return "CRAM: cannot open " + path;

    u32 hdr[2] = {0, 0};
    if (std::fread(hdr, sizeof(u32), 2, f) != 2 || hdr[0] != MAGIC || hdr[1] != u32(NPARAM))
    {
        std::fclose(f);
        return "CRAM: not a compatible net (bad magic or parameter count, expected "
             + std::to_string(NPARAM) + ")";
    }

    std::vector<float> w(NPARAM);
    const size_t       got = std::fread(w.data(), sizeof(float), NPARAM, f);
    unsigned char      extra;
    const bool         trailing = std::fread(&extra, 1, 1, f) == 1;
    std::fclose(f);

    if (got != size_t(NPARAM) || trailing)
        return "CRAM: wrong file size";
    for (float v : w)
        if (!std::isfinite(v))
            return "CRAM: net contains NaN/Inf";

    g_w      = std::move(w);
    g_loaded = true;
    return "CRAM net loaded: " + path;
}

void add_ordering_bonus(const Position& pos, ExtMove* first, ExtMove* last) {
    const int n = int(last - first);
    if (!g_loaded || n <= 1 || g_weight == 0)
        return;

    Ctx   ctx(pos);
    float sc[MAX_MOVES], mean = 0.f;
    for (int i = 0; i < n; ++i)
        mean += (sc[i] = ctx.score(first[i]));
    mean /= float(n);

    for (int i = 0; i < n; ++i)
        first[i].value += std::clamp(int(g_weight * (sc[i] - mean)), -32768, 32768);
}

bool gate_pick(const Position& pos, const Move* cands, int n, Move& best, float& confidence) {
    if (!g_loaded || n <= 0)
        return false;

    Ctx   ctx(pos);
    float sc[MAX_MOVES];
    int   b = 0;
    for (int i = 0; i < n; ++i)
    {
        sc[i] = ctx.score(cands[i]);
        if (sc[i] > sc[b])
            b = i;
    }
    float z = 0.f;
    for (int i = 0; i < n; ++i)
        z += std::exp(sc[i] - sc[b]);

    best       = cands[b];
    confidence = 1.0f / z;
    return true;
}

std::string report(const Position& pos, bool chess960) {
    if (!g_loaded)
        return "cram: no net loaded (setoption name CRAMFile value <path>)";

    Ctx                                ctx(pos);
    std::vector<std::pair<float, Move>> v;
    for (const auto& m : MoveList<ALL>(pos))
        v.emplace_back(ctx.score(m), m);
    if (v.empty())
        return "cram: no legal moves";
    std::sort(v.begin(), v.end(), [](auto& a, auto& b) { return a.first > b.first; });

    float z = 0.f;
    for (auto& e : v)
        z += std::exp(e.first - v[0].first);

    std::ostringstream os;
    os << std::fixed << std::setprecision(6);
    for (auto& e : v)
        os << "cram_move " << UCIEngine::move(e.second, chess960) << " " << e.first << "\n";
    os << "cram_confidence " << 1.0f / z << "\n";
    return os.str();
}

}  // namespace Stockfish::CRAM
