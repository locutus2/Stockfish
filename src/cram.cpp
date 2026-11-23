/*
  CRAM move prior for Stockfish (experimental add-on, not part of upstream Stockfish)

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/

#include "cram.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>

#include "bitboard.h"
#include "misc.h"
#include "movegen.h"
#include "position.h"
#include "uci.h"

// The default net is compiled into the binary as a plain C++ array, so that no assembler
// tricks (incbin) or external files are needed and the engine behaves the same on every
// platform, compiler and build type (also on Fishtest workers, which only get the binary).
// cram_net.inc is generated from model-i8.bin by scripts/gen_cram_inc.py; regenerate it
// whenever model-i8.bin changes.
namespace {
const unsigned char gEmbeddedCRAMData[] = {
#include "cram_net.inc"
};
const unsigned int gEmbeddedCRAMSize = sizeof(gEmbeddedCRAMData);
}  // namespace

namespace Stockfish::CRAM {

bool g_loaded = false;
int  g_plies  = 256;

namespace {

std::vector<int8_t> g_w;  // the net, NPARAM int8 weights; only replaced while no search runs
i32                 g_mult[NGROUPS_H]; // value of weight i of group g = g_w[i] * g_mult[g] / 2^FRAC
bool                g_hidden = false;  // net with hidden layer (CRM4)?
i64                 g_inv[NZ][64];     // g_inv[k][d] = 2^(ZFRAC+FRAC) / (mass_k + d * 2^FRAC)
int                g_weight = 8192;
int                g_gate   = 0;

// Fixed-point exp: e^(-d / 2^FRAC) with 30 fractional bits (1.0 = 2^30), integer arithmetic only (d >= 0).
// Table 2^(-k/256) * 2^30, k = 0..256, with linear integer interpolation; max error ~1.4e-6.
// Must stay identical to chess_cram.cpp.
constexpr i32 EXP2_TAB[257] = {
    1073741824, 1070838486, 1067942999, 1065055341, 1062175491, 1059303428, 1056439131, 1053582579,
    1050733751, 1047892626, 1045059183, 1042233401, 1039415261, 1036604740, 1033801819, 1031006477,
    1028218693, 1025438448, 1022665720, 1019900489, 1017142735, 1014392438, 1011649578, 1008914134,
    1006186087, 1003465416, 1000752102, 998046124, 995347464, 992656100, 989972014, 987295185,
    984625594, 981963222, 979308048, 976660054, 974019220, 971385527, 968758955, 966139485,
    963527098, 960921775, 958323496, 955732243, 953147997, 950570738, 948000448, 945437108,
    942880699, 940331203, 937788600, 935252872, 932724001, 930201967, 927686753, 925178340,
    922676710, 920181844, 917693724, 915212331, 912737649, 910269657, 907808339, 905353676,
    902905651, 900464244, 898029440, 895601218, 893179563, 890764456, 888355878, 885953814,
    883558244, 881169153, 878786521, 876410331, 874040567, 871677210, 869320244, 866969651,
    864625413, 862287515, 859955938, 857630665, 855311680, 852998965, 850692504, 848392279,
    846098274, 843810471, 841528855, 839253408, 836984114, 834720956, 832463917, 830212982,
    827968132, 825729353, 823496627, 821269938, 819049271, 816834607, 814625932, 812423229,
    810226483, 808035676, 805850792, 803671817, 801498734, 799331526, 797170178, 795014675,
    792865000, 790721137, 788583072, 786450787, 784324269, 782203500, 780088465, 777979150,
    775875538, 773777614, 771685363, 769598769, 767517817, 765442492, 763372778, 761308661,
    759250125, 757197155, 755149737, 753107854, 751071493, 749040637, 747015274, 744995386,
    742980960, 740971982, 738968435, 736970306, 734977579, 732990241, 731008277, 729031671,
    727060411, 725094480, 723133865, 721178552, 719228525, 717283772, 715344277, 713410026,
    711481005, 709557200, 707638598, 705725183, 703816941, 701913860, 700015924, 698123120,
    696235434, 694352853, 692475362, 690602947, 688735596, 686873293, 685016026, 683163781,
    681316545, 679474303, 677637043, 675804750, 673977412, 672155015, 670337545, 668524990,
    666717336, 664914570, 663116678, 661323648, 659535466, 657752119, 655973594, 654199878,
    652430958, 650666822, 648907455, 647152846, 645402981, 643657847, 641917433, 640181724,
    638450708, 636724373, 635002706, 633285695, 631573326, 629865587, 628162466, 626463950,
    624770026, 623080683, 621395908, 619715688, 618040012, 616368866, 614702239, 613040119,
    611382493, 609729349, 608080675, 606436459, 604796689, 603161352, 601530438, 599903933,
    598281827, 596664106, 595050760, 593441776, 591837143, 590236848, 588640881, 587049229,
    585461881, 583878825, 582300049, 580725543, 579155293, 577589290, 576027521, 574469975,
    572916640, 571367506, 569822560, 568281792, 566745190, 565212742, 563684439, 562160268,
    560640218, 559124278, 557612438, 556104685, 554601009, 553101399, 551605844, 550114332,
    548626854, 547143398, 545663953, 544188508, 542717053, 541249576, 539786068, 538326517,
    536870912,
};
inline i64 expNeg(i64 d) {
    if (d >= (i64(32) << FRAC)) return 0;                    // e^-32 ~ 1e-14: below the resolution
    const i64 y  = (d * 1512775) >> FRAC;                    // d * log2(e); 1512775 = round(log2(e) * 2^20)
    const int     n  = int(y >> FRAC);                           // integer part of the exponent
    const i64 fr = y & ((i64(1) << FRAC) - 1);
    const int     ix = int(fr >> 12);
    const i64 rm = fr & 4095;
    const i64 v  = (i64(EXP2_TAB[ix]) * (4096 - rm) + i64(EXP2_TAB[ix + 1]) * rm) >> 12;
    return v >> n;
}

// Direction class of a move between oriented squares (must match dirClass() in chess_cram.cpp):
// 0..7 queen directions, 8 knight jump, 9 stay.
struct DirTable {
    uint8_t t[64][64];
    DirTable() {
        for (int of = 0; of < 64; ++of)
            for (int ot = 0; ot < 64; ++ot)
            {
                const int df = (ot & 7) - (of & 7), dr = (ot >> 3) - (of >> 3);
                int       c;
                if (of == ot)
                    c = 9;
                else if ((std::abs(df) == 1 && std::abs(dr) == 2)
                         || (std::abs(df) == 2 && std::abs(dr) == 1))
                    c = 8;
                else
                {
                    c = ((df > 0) - (df < 0) + 1) * 3 + ((dr > 0) - (dr < 0) + 1);
                    c = c < 4 ? c : c - 1;
                }
                t[of][ot] = uint8_t(c);
            }
    }
};
const DirTable g_dir;

// Per-position evaluation context. Field scores are computed lazily and shared by all moves
// that start/end on the same square (about 2x fewer dot products than moves).
// All scores are i64 fixed-point values in units of 2^-FRAC.
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

        if (g_hidden)
            hidden_setup();
    }

    // ---- hidden layer (see chess_cram.cpp, hidInt(), for the specification) ----------------------
    // Works in real squares for the relaxation (the graph is the same up to relabeling the squares) and
    // in side-to-move oriented squares (square ^ flip) for the direction classes.
    void hidden_setup() {
        const Bitboard occ = pos.pieces();
        for (Bitboard b = occ; b;)
        {
            const Square s  = pop_lsb(b);
            const Piece  pc = pos.piece_on(s);
            const Bitboard a = type_of(pc) == PAWN ? Attacks::attacks_bb(PAWN, s, color_of(pc))
                                                   : Attacks::attacks_bb(type_of(pc), s, occ);
            adj[s] |= a;
            for (Bitboard t = a; t;)
                adj[pop_lsb(t)] |= square_bb(s);
        }
        const int goal[NZ] = {int(pos.square<KING>(~us)), int(pos.square<KING>(us))};

        i64 za[NZ][64] = {}, zb[NZ][64];
        i64(*cur)[64]  = za;
        i64(*nxt)[64]  = zb;
        for (int t = 0; t < HT; ++t)
        {
            for (int i = 0; i < 64; ++i)
            {
                const int d  = popcount(adj[i]);
                i64       s0 = i == goal[0] ? (i64(1) << ZFRAC) : 0;
                i64       s1 = i == goal[1] ? (i64(1) << ZFRAC) : 0;
                for (Bitboard b = adj[i]; b;)
                {
                    const Square j = pop_lsb(b);
                    s0 += cur[0][j];
                    s1 += cur[1][j];
                }
                nxt[0][i] = (s0 * g_inv[0][d]) >> ZFRAC;
                nxt[1][i] = (s1 * g_inv[1][d]) >> ZFRAC;
            }
            std::swap(cur, nxt);
        }
        for (int k = 0; k < NZ; ++k)
        {
            i64 zmax = 0;
            for (int i = 0; i < 64; ++i)
            {
                z[k][i] = cur[k][i];
                zmax    = std::max(zmax, z[k][i]);
            }
            zrecip[k] = zmax > 0 ? (i64(1) << 62) / zmax : 0;
        }
    }

    // zhat = z / max z with 16 fraction bits: (z * rz) >> 46 with rz = 2^62 / max z (one division per field).
    i64 zhat(int k, int sq) { return (z[k][sq] * zrecip[k]) >> 46; }

    i64 attention_logit(int k, int dir, int sq) {
        return (i64(g_mult[G_BETA]) * g_w[OFF_BETA + k * NDIR + dir] * zhat(k, sq)) / 65536;
    }

    // Attention window of the piece on oriented square 'of': "stay", every attacked square, the pawn
    // pushes and the castling target squares. Computes max logit and 1 / sum(exp) per field.
    void window(int of) {
        const Square from = Square(of ^ flip);
        const Piece  pc   = pos.piece_on(from);
        const PieceType pt = type_of(pc);
        Bitboard bb = pt == PAWN ? Attacks::attacks_bb(PAWN, from, us)
                                 : Attacks::attacks_bb(pt, from, pos.pieces());
        bb |= square_bb(from);
        if (pt == PAWN)
        {
            const Square push = from + (us == WHITE ? NORTH : SOUTH);
            bb |= square_bb(push);
            if (relative_rank(us, from) == RANK_2)
                bb |= square_bb(push + (us == WHITE ? NORTH : SOUTH));
        }
        else if (pt == KING && relative_rank(us, from) == RANK_1)
            bb |= square_bb(make_square(FILE_G, rank_of(from))) | square_bb(make_square(FILE_C, rank_of(from)));

        i64 a[NZ][32];
        int n = 0;
        for (Bitboard b = bb; b; ++n)
        {
            const Square s   = pop_lsb(b);
            const int    dir = g_dir.t[of][int(s) ^ flip];
            for (int k = 0; k < NZ; ++k)
                a[k][n] = attention_logit(k, dir, s);
        }
        for (int k = 0; k < NZ; ++k)
        {
            i64 amax = a[k][0], wsum = 0;
            for (int u = 1; u < n; ++u)
                amax = std::max(amax, a[k][u]);
            for (int u = 0; u < n; ++u)
                wsum += expNeg(amax - a[k][u]);
            winAmax[k][of] = amax;
            winRinv[k][of] = (i64(1) << 60) / wsum;
        }
        winDone |= u64(1) << of;
    }

    // Contribution of the hidden layer to the score of a move of piece type 'type' from oriented square
    // 'of' to the real square 'to' (units 2^-FRAC).
    i64 hidden_score(int of, Square to, int type) {
        if (!(winDone >> of & 1))
            window(of);
        const int dir = g_dir.t[of][int(to) ^ flip];
        i64       r   = 0;
        for (int k = 0; k < NZ; ++k)
        {
            const i64 w     = expNeg(winAmax[k][of] - attention_logit(k, dir, to));
            const i64 alpha = (w * winRinv[k][of]) >> 30;
            r += (i64(g_mult[G_KAPPA]) * g_w[OFF_KAPPA + k * 6 + type] * alpha) / (i64(1) << 30);
        }
        return r;
    }


    // Field score = bias + sum of the piece-square weights: int8 weights are summed exactly in
    // int32 (at most 32 pieces), then each group's integer multiplier is applied once.
    i64 field(int base, int bias, int gW, int gB, int o) const {
        const int8_t* w   = g_w.data();
        const int8_t* row = w + base + o * NF;
        int           acc = 0;
        for (int j = 0; j < nfeat; ++j)
            acc += row[fidx[j]];
        return i64(g_mult[gW]) * acc + i64(g_mult[gB]) * w[bias + o];
    }

    i64 score(Move m) {
        const int8_t* w    = g_w.data();
        Square        from = m.from_sq(), to = m.to_sq();
        const bool    castling = m.type_of() == CASTLING;
        if (castling)  // Stockfish encodes castling as "king takes own rook"; use the king's target
            to = make_square(to > from ? FILE_G : FILE_C, rank_of(from));

        const int of = int(from) ^ flip, ot = int(to) ^ flip;
        if (!(doneF >> of & 1))
            fF[of] = field(OFF_WF, OFF_BF, G_WF, G_BF, of), doneF |= u64(1) << of;
        if (!(doneT >> ot & 1))
            fT[ot] = field(OFF_WT, OFF_BT, G_WT, G_BT, ot), doneT |= u64(1) << ot;

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

        return (g_hidden ? hidden_score(of, to, type) : 0) + fF[of] + fT[ot]
             + i64(g_mult[G_K]) * w[OFF_K + of * 64 + ot]
             + (promo ? i64(g_mult[G_P]) * w[OFF_P + promo] : 0)
             + i64(g_mult[G_KMC]) * w[OFF_KMC + type * 7 + cap]
             + i64(g_mult[G_A]) * w[OFF_A + (type * 4 + eN) * 4 + oN]
             + i64(g_mult[G_AF]) * w[OFF_AF + type * 4 + eF];
    }

    const Position& pos;
    Color           us;
    int             flip;
    int             nfeat = 0;
    u16             fidx[32];
    u64             doneF = 0, doneT = 0;
    i64             fF[64], fT[64];
    // hidden layer state
    Bitboard adj[64] = {};         // piece-mobility graph (real squares)
    i64      z[NZ][64], zrecip[NZ];
    u64      winDone = 0;
    i64      winAmax[NZ][64], winRinv[NZ][64];
};

// Softmax probability of the best move in millionths, integers only.
i64 confidence_ppm(const i64* sc, int n, int best) {
    i64 sum = 0;
    for (int i = 0; i < n; ++i)
        sum += expNeg(sc[best] - sc[i]);
    return (i64(1000000) << 30) / sum;
}

}  // namespace

bool loaded() { return g_loaded; }
bool has_hidden() { return g_loaded && g_hidden; }

void unload() {
    g_loaded = false;
    g_hidden = false;
    std::vector<int8_t>().swap(g_w);
}

void set_weight(int w) { g_weight = w; }
void set_plies(int p) { g_plies = p; }
void set_gate(int pct) { g_gate = pct; }
int  gate_percent() { return g_gate; }

namespace {

// Parses a CRM3 or CRM4 net from memory and installs it. Returns an error message, or an empty string.
std::string parse_net(const unsigned char* data, size_t size) {
    if (size < 2 * sizeof(u32))
        return "cannot read header";

    u32 magic = 0;
    std::memcpy(&magic, data, sizeof(u32));  // little-endian host assumed, like the NNUE loader
    if (magic == MAGIC_F32 || magic == MAGIC_Q8F)
        return "old net with float weights/scales (CRM1/CRM2) is no longer supported; "
               "retrain with the integer chess_cram";

    const bool hid = magic == MAGIC_H;
    if (magic != MAGIC && !hid)
        return "not a compatible int8 net (bad magic; expected CRM3 or CRM4)";

    const int    np      = hid ? NPARAM_H : NPARAM;
    const int    ng      = hid ? NGROUPS_H : NGROUPS;
    const size_t nhdr    = hid ? 6 : 4;
    const size_t HEADER  = nhdr * sizeof(u32) + ng * sizeof(i32) + (hid ? NZ * sizeof(i32) : 0);
    if (size < HEADER)
        return "wrong file size";

    u32 hdr[6] = {};
    std::memcpy(hdr, data, nhdr * sizeof(u32));
    if (hdr[1] != u32(np) || hdr[2] != u32(ng) || hdr[3] != u32(FRAC)
        || (hid && (hdr[4] != u32(NZ) || hdr[5] != u32(HT))))
        return "not a compatible int8 net (bad parameter count, group count, fixed-point format or "
               "hidden-layer shape; expected "
             + std::to_string(np) + " parameters)";
    if (size != HEADER + size_t(np))
        return "wrong file size";

    i32 mult[NGROUPS_H] = {}, mass[NZ] = {};
    std::memcpy(mult, data + nhdr * sizeof(u32), ng * sizeof(i32));
    for (int g = 0; g < ng; ++g)
        if (mult[g] < 1 || mult[g] > MAX_MULT)
            return "invalid multiplier";
    if (hid)
    {
        std::memcpy(mass, data + nhdr * sizeof(u32) + ng * sizeof(i32), sizeof(mass));
        for (i32 v : mass)
            if (v < MASS_MIN_Q || v > MAX_MASS_Q)
                return "invalid hidden-layer mass";
    }

    const int8_t* w = reinterpret_cast<const int8_t*>(data + HEADER);
    g_w.assign(w, w + np);
    std::fill(g_mult, g_mult + NGROUPS_H, 0);
    std::copy(mult, mult + ng, g_mult);
    // A hidden layer whose kappa weights are all zero adds exactly 0 to every score: skip it (same
    // result, no cost). This also makes a CRM3 net converted with scripts/crm3_to_crm4.py behave
    // exactly like the original.
    g_hidden = hid && std::any_of(w + OFF_KAPPA, w + OFF_BETA, [](int8_t v) { return v != 0; });
    if (g_hidden)
        for (int k = 0; k < NZ; ++k)
            for (int d = 0; d < 64; ++d)
                g_inv[k][d] = (i64(1) << (ZFRAC + FRAC)) / (i64(mass[k]) + (i64(d) << FRAC));
    g_loaded = true;
    return "";
}

}  // namespace

// Loads a net. The default file name selects the net embedded in the binary (like the NNUE
// net), so that the engine works without any external file, e.g. on Fishtest workers.
// Any other name is read from disk.
std::optional<std::string> load(const std::string& path) {
    if (path.empty() || path == "<empty>")
    {
        unload();
        return "CRAM net unloaded";
    }

    if (path == DefaultFile && gEmbeddedCRAMSize > 1)
    {
        const std::string err = parse_net(gEmbeddedCRAMData, gEmbeddedCRAMSize);
        if (!err.empty())
            return "CRAM: embedded net is invalid: " + err;
        return "CRAM int8 net loaded: embedded " + std::string(DefaultFile) + (g_hidden ? " (with hidden layer)" : "");
    }

    // The option value is UTF-8; path_from_utf8 makes this work with non-ASCII paths on Windows.
    std::ifstream f(path_from_utf8(path), std::ios::binary);
    if (!f)
        return "CRAM: cannot open " + path;

    std::vector<unsigned char> buf;
    char                       chunk[65536];
    while (f.read(chunk, sizeof(chunk)) || f.gcount() > 0)
    {
        buf.insert(buf.end(), chunk, chunk + f.gcount());
        if (buf.size() > (size_t(1) << 26))  // far larger than any valid net
            break;
    }

    const std::string err = parse_net(buf.data(), buf.size());
    if (!err.empty())
        return "CRAM: " + err;
    return "CRAM int8 net loaded: " + path + (g_hidden ? " (with hidden layer)" : "");
}

void add_ordering_bonus(const Position& pos, ExtMove* first, ExtMove* last, int /*depth*/) {
    const int n = int(last - first);
    if (!g_loaded || n <= 1 || g_weight == 0)
        return;

    Ctx ctx(pos);
    i64 sc[MAX_MOVES], sum = 0;
    for (int i = 0; i < n; ++i)
        sum += (sc[i] = ctx.score(first[i]));
    const i64 mean = sum / n;  // integer division truncates toward zero: well defined

    for (int i = 0; i < n; ++i)
    {
        // weight * (score - mean), converted from 2^-FRAC fixed point; truncation toward zero
        const i64 bonus = i64(g_weight) * (sc[i] - mean) / (i64(1) << FRAC);
        first[i].value += int(std::clamp<i64>(bonus, -32768, 32768));
        //first[i].value += int(std::clamp<i64>(bonus / std::max(depth, 1), -32768, 32768));
    }
}

bool gate_pick(const Position& pos, const Move* cands, int n, Move& best, i64& confidencePpm) {
    if (!g_loaded || n <= 0)
        return false;

    Ctx ctx(pos);
    i64 sc[MAX_MOVES];
    int b = 0;
    for (int i = 0; i < n; ++i)
    {
        sc[i] = ctx.score(cands[i]);
        if (sc[i] > sc[b] || (sc[i] == sc[b] && cands[i].raw() < cands[b].raw()))
            b = i;
    }

    best          = cands[b];
    confidencePpm = confidence_ppm(sc, n, b);
    return true;
}

std::string report(const Position& pos, bool chess960) {
    if (!g_loaded)
        return "cram: no net loaded (setoption name CRAMFile value <path>)";

    Ctx                              ctx(pos);
    std::vector<std::pair<i64, Move>> v;
    for (const auto& m : MoveList<ALL>(pos))
        v.emplace_back(ctx.score(m), m);
    if (v.empty())
        return "cram: no legal moves";
    // Total order: score descending, then Move::raw() ascending
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
        return a.first != b.first ? a.first > b.first : a.second.raw() < b.second.raw();
    });

    std::vector<i64> sc;
    for (auto& e : v)
        sc.push_back(e.first);

    // Scores are integers in units of 2^-FRAC; confidence in millionths.
    std::ostringstream os;
    for (auto& e : v)
        os << "cram_move " << UCIEngine::move(e.second, chess960) << " " << e.first << "\n";
    os << "cram_confidence_ppm " << confidence_ppm(sc.data(), int(sc.size()), 0) << "\n";
    return os.str();
}

}  // namespace Stockfish::CRAM
