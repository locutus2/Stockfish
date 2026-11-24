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

#if defined(USE_AVX2)
    #include <immintrin.h>
#elif defined(USE_SSE41)
    #include <smmintrin.h>
#elif defined(USE_NEON)
    #include <arm_neon.h>
#endif

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

// SIMD is used for the piece-square sums where the build supports it. All arithmetic is integer
// and exact, so every code path (SIMD or scalar) produces bit-identical scores.
#if defined(USE_AVX2) || defined(USE_SSE41) || defined(USE_NEON)
    #define CRAM_SIMD 1
#else
    #define CRAM_SIMD 0
#endif

namespace {

std::vector<int8_t> g_w;  // the net, NPARAM int8 weights; only replaced while no search runs
#if CRAM_SIMD
// The piece-square weights transposed to [feature][square] (WF block, then WT block), so that
// adding one active feature updates the sums of all 64 squares with contiguous loads.
std::vector<int8_t> g_wt;
#endif
i32                 g_mult[NGROUPS_H]; // value of weight i of group g = g_w[i] * g_mult[g] / 2^FRAC
bool                g_hidden = false;  // net with hidden layer (CRM4)?
i64                 g_inv[NZ][64];     // g_inv[k][d] = 2^(ZFRAC+FRAC) / (mass_k + d * 2^FRAC)
i64                 g_bq[NZ][NDIR];    // beta * its scale:  g_mult[G_BETA] * q_beta[field][direction class]
i32                 g_bq32[NZ][16];    // g_bq as int32 (|value| <= 127 * 2^24 < 2^31), padded for gathers
i64                 g_kc[NZ][6];       // kappa * its scale: g_mult[G_KAPPA] * q_kappa[field][piece type]
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

#if defined(USE_AVX2)
// g_expPair.v[i] = EXP2_TAB[i] | EXP2_TAB[i + 1] << 32: lets the vector code fetch both table entries
// that expNeg() interpolates with one 64-bit gather.
struct ExpPair {
    alignas(32) u64 v[256];
    ExpPair() {
        for (int i = 0; i < 256; ++i)
            v[i] = u64(u32(EXP2_TAB[i])) | u64(u32(EXP2_TAB[i + 1])) << 32;
    }
};
const ExpPair g_expPair;
#endif

// floor(2^60 / d) for 2^30 <= d < 2^40, exactly, without a 64-bit integer division (which takes 40-90 cycles
// on many CPUs): a double division is correct to within 1 of the result, an integer check fixes the rest.
inline i64 recip_2_60(i64 d) {
    i64 q = i64(1152921504606846976.0 / double(d));
    i64 r = (i64(1) << 60) - q * d;  // remainder of the estimate
    if (r < 0)
        --q;
    else if (r >= d)
        ++q;
    return q;
}

// Number of attackers (saturated at 3) of every square, kept as three bit planes.
struct AttackCount {
    Bitboard c1 = 0, c2 = 0, c3 = 0;  // squares attacked at least once / twice / three times

    void add(Bitboard attacks) {
        c3 |= c2 & attacks;
        c2 |= c1 & attacks;
        c1 |= attacks;
    }
    int at(Square s) const { return int(c1 >> s & 1) + int(c2 >> s & 1) + int(c3 >> s & 1); }
};

// With fewer moves than this, summing only the needed squares is cheaper than summing all 64.
constexpr int SIMD_MIN_MOVES = 8;

#if CRAM_SIMD
// sumF[o] / sumT[o] = sum of the int8 weights of all active features for from-square / to-square
// o. At most 32 features of magnitude <= 127 are active, so the sums (<= 4064) fit into int16.
inline void piece_square_sums(const u16* fidx, int nfeat, int16_t* sumF, int16_t* sumT) {
    const int8_t* wf = g_wt.data();
    const int8_t* wt = wf + size_t(64) * NF;
    #if defined(USE_AVX2)
    __m256i f0 = _mm256_setzero_si256(), f1 = f0, f2 = f0, f3 = f0;
    __m256i t0 = f0, t1 = f0, t2 = f0, t3 = f0;
    for (int j = 0; j < nfeat; ++j)
    {
        const __m128i* rf = reinterpret_cast<const __m128i*>(wf + fidx[j] * 64);
        const __m128i* rt = reinterpret_cast<const __m128i*>(wt + fidx[j] * 64);
        f0 = _mm256_add_epi16(f0, _mm256_cvtepi8_epi16(_mm_loadu_si128(rf + 0)));
        f1 = _mm256_add_epi16(f1, _mm256_cvtepi8_epi16(_mm_loadu_si128(rf + 1)));
        f2 = _mm256_add_epi16(f2, _mm256_cvtepi8_epi16(_mm_loadu_si128(rf + 2)));
        f3 = _mm256_add_epi16(f3, _mm256_cvtepi8_epi16(_mm_loadu_si128(rf + 3)));
        t0 = _mm256_add_epi16(t0, _mm256_cvtepi8_epi16(_mm_loadu_si128(rt + 0)));
        t1 = _mm256_add_epi16(t1, _mm256_cvtepi8_epi16(_mm_loadu_si128(rt + 1)));
        t2 = _mm256_add_epi16(t2, _mm256_cvtepi8_epi16(_mm_loadu_si128(rt + 2)));
        t3 = _mm256_add_epi16(t3, _mm256_cvtepi8_epi16(_mm_loadu_si128(rt + 3)));
    }
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(sumF) + 0, f0);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(sumF) + 1, f1);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(sumF) + 2, f2);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(sumF) + 3, f3);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(sumT) + 0, t0);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(sumT) + 1, t1);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(sumT) + 2, t2);
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(sumT) + 3, t3);
    #elif defined(USE_SSE41)
    // 64 int16 sums = 8 SSE registers per table; the two tables are processed one after another.
    for (int table = 0; table < 2; ++table)
    {
        const int8_t* w   = table ? wt : wf;
        int16_t*      out = table ? sumT : sumF;
        __m128i       a0 = _mm_setzero_si128(), a1 = a0, a2 = a0, a3 = a0, a4 = a0, a5 = a0,
                      a6 = a0, a7 = a0;
        for (int j = 0; j < nfeat; ++j)
        {
            const __m128i* r = reinterpret_cast<const __m128i*>(w + fidx[j] * 64);
            // _mm_cvtepi8_epi16 sign-extends the low 8 bytes of its argument
            a0 = _mm_add_epi16(a0, _mm_cvtepi8_epi16(_mm_loadu_si128(r + 0)));
            a1 = _mm_add_epi16(a1, _mm_cvtepi8_epi16(_mm_srli_si128(_mm_loadu_si128(r + 0), 8)));
            a2 = _mm_add_epi16(a2, _mm_cvtepi8_epi16(_mm_loadu_si128(r + 1)));
            a3 = _mm_add_epi16(a3, _mm_cvtepi8_epi16(_mm_srli_si128(_mm_loadu_si128(r + 1), 8)));
            a4 = _mm_add_epi16(a4, _mm_cvtepi8_epi16(_mm_loadu_si128(r + 2)));
            a5 = _mm_add_epi16(a5, _mm_cvtepi8_epi16(_mm_srli_si128(_mm_loadu_si128(r + 2), 8)));
            a6 = _mm_add_epi16(a6, _mm_cvtepi8_epi16(_mm_loadu_si128(r + 3)));
            a7 = _mm_add_epi16(a7, _mm_cvtepi8_epi16(_mm_srli_si128(_mm_loadu_si128(r + 3), 8)));
        }
        __m128i* o = reinterpret_cast<__m128i*>(out);
        _mm_storeu_si128(o + 0, a0);
        _mm_storeu_si128(o + 1, a1);
        _mm_storeu_si128(o + 2, a2);
        _mm_storeu_si128(o + 3, a3);
        _mm_storeu_si128(o + 4, a4);
        _mm_storeu_si128(o + 5, a5);
        _mm_storeu_si128(o + 6, a6);
        _mm_storeu_si128(o + 7, a7);
    }
    #elif defined(USE_NEON)
    // 64 int16 sums = 8 NEON registers per table; vaddw_s8 adds a sign-extended int8 vector.
    int16x8_t f[8], t[8];
    for (int k = 0; k < 8; ++k)
        f[k] = t[k] = vdupq_n_s16(0);
    for (int j = 0; j < nfeat; ++j)
    {
        const int8_t* rf = wf + fidx[j] * 64;
        const int8_t* rt = wt + fidx[j] * 64;
        for (int q = 0; q < 4; ++q)
        {
            const int8x16_t vf = vld1q_s8(rf + 16 * q);
            const int8x16_t vt = vld1q_s8(rt + 16 * q);
            f[2 * q]           = vaddw_s8(f[2 * q], vget_low_s8(vf));
            f[2 * q + 1]       = vaddw_s8(f[2 * q + 1], vget_high_s8(vf));
            t[2 * q]           = vaddw_s8(t[2 * q], vget_low_s8(vt));
            t[2 * q + 1]       = vaddw_s8(t[2 * q + 1], vget_high_s8(vt));
        }
    }
    for (int k = 0; k < 8; ++k)
    {
        vst1q_s16(sumF + 8 * k, f[k]);
        vst1q_s16(sumT + 8 * k, t[k]);
    }
    #endif
}
#endif  // CRAM_SIMD

// Direction class of a move between oriented squares (must match dirClass() in chess_cram.cpp):
// 0..7 queen directions, 8 knight jump, 9 stay.
struct DirTable {
    uint8_t t[64][64];       // [oriented from][oriented to]
    uint8_t row[2][64][64];  // [flip != 0][real from][real to]: the same class, indexed by real squares
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
        for (int f = 0; f < 2; ++f)
            for (int from = 0; from < 64; ++from)
                for (int to = 0; to < 64; ++to)
                    row[f][from][to] = t[from ^ (f ? 56 : 0)][to ^ (f ? 56 : 0)];
    }
};
const DirTable g_dir;

// Per-position evaluation context. Field scores are computed lazily and shared by all moves
// that start/end on the same square (about 2x fewer dot products than moves). With enough moves
// the piece-square sums of all squares are instead computed at once with SIMD.
// All scores are i64 fixed-point values in units of 2^-FRAC.
struct Ctx {
    Ctx(const Position& p, int numMoves) :
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

#if CRAM_SIMD
        if (numMoves >= SIMD_MIN_MOVES)
        {
            piece_square_sums(fidx, nfeat, sumF, sumT);
            haveSums = true;
        }
#else
        (void) numMoves;
#endif

        // Attackers of every square (saturated at 3) for the opponent and for us. This equals
        // popcount(pos.attackers_to(s) & pieces(color)) clamped to 3, but is computed once per
        // position from the attack sets of all pieces instead of twice per move.
        const Bitboard occ = pos.pieces();
        for (Color c : {WHITE, BLACK})
        {
            AttackCount&   a     = c == us ? ours : theirs;
            const Bitboard pawns = pos.pieces(c, PAWN);
            a.add(shift(pawns, c == WHITE ? NORTH_WEST : SOUTH_WEST));
            a.add(shift(pawns, c == WHITE ? NORTH_EAST : SOUTH_EAST));
            if (g_hidden)  // the hidden layer needs the attack set of every single pawn
                for (Bitboard pw = pawns; pw;)
                {
                    const Square s = pop_lsb(pw);
                    add_edges(s, Attacks::attacks_bb(PAWN, s, c));
                }
            for (PieceType pt : {KNIGHT, BISHOP, ROOK, QUEEN})
            {
                Bitboard pcs = pos.pieces(c, pt);
                while (pcs)
                {
                    const Square   s  = pop_lsb(pcs);
                    const Bitboard at = Attacks::attacks_bb(pt, s, occ);
                    a.add(at);
                    if (g_hidden)
                        add_edges(s, at);
                }
            }
            const Square   ksq = pos.square<KING>(c);
            const Bitboard kat = Attacks::attacks_bb(KING, ksq, occ);
            a.add(kat);
            if (g_hidden)
                add_edges(ksq, kat);
        }

        if (g_hidden)
            hidden_setup();
    }

    // ---- hidden layer (see chess_cram.cpp, hidInt(), for the specification) ----------------------
    // The code below is an optimized but bit-identical implementation of that specification:
    //  * The relaxation works in real squares (relabeling the squares does not change the graph) and
    //    exploits that z = 0 at the start: step 1 only touches the goal square, step 2 its neighbours,
    //    and so on. Every step scatters only from the squares that are non-zero and finishes only the
    //    squares that received something, instead of looping over all 64 squares.
    //  * The attention weights of a window are computed once per from-square (exp of all entries at
    //    once, with SIMD where available) and stored per destination, so a move only needs a lookup.
    // Adds the edges of a piece on s that attacks the squares 'at' (the graph is undirected).
    void add_edges(Square s, Bitboard at) {
        adj[s] |= at;
        for (Bitboard t = at; t;)
            adj[pop_lsb(t)] |= square_bb(s);
    }

    void hidden_setup() {
        std::memset(slotMap, -1, sizeof(slotMap));
        std::memset(zhv, 0, sizeof(zhv));

        for (int k = 0; k < NZ; ++k)
        {
            const Square goal = pos.square<KING>(k == 0 ? ~us : us);  // field 0: opponent king, 1: own king
            i64          zc[64];  // only read for squares in 'active', which have been written before
            i64          sum[64];
            Bitboard     active = 0;  // squares that may be non-zero
            for (int t = 0; t < HT; ++t)
            {
                Bitboard touched = square_bb(goal);
                for (Bitboard b = active; b;)
                    touched |= adj[pop_lsb(b)];
                for (Bitboard b = touched; b;)
                    sum[pop_lsb(b)] = 0;
                sum[goal] = i64(1) << ZFRAC;
                for (Bitboard b = active; b;)
                {
                    const Square j = pop_lsb(b);
                    const i64    v = zc[j];
                    for (Bitboard n = adj[j]; n;)
                        sum[pop_lsb(n)] += v;
                }
                for (Bitboard b = active; b;)  // squares without a non-zero neighbour become 0
                    zc[pop_lsb(b)] = 0;
                for (Bitboard b = touched; b;)
                {
                    const Square i = pop_lsb(b);
                    zc[i]          = (sum[i] * g_inv[k][popcount(adj[i])]) >> ZFRAC;
                }
                active = touched;
            }
            i64 zmax = 0;
            for (Bitboard b = active; b;)
                zmax = std::max(zmax, zc[pop_lsb(b)]);
            const i64 zrecip = zmax > 0 ? (i64(1) << 62) / zmax : 0;
            // zhat = z / max z with 16 fraction bits: (z * zrecip) >> 46, zrecip = 2^62 / max z.
            // Squares outside 'active' have z = 0 (zhv is zero-initialized).
            for (Bitboard b = active; b;)
            {
                const Square i = pop_lsb(b);
                zhv[k][i]      = i32((zc[i] * zrecip) >> 46);
            }
        }
    }

    // Squares of the attention window of the piece on 'from': "stay", every attacked square, the pawn
    // pushes and the castling target squares.
    Bitboard window_squares(Square from) const {
        const PieceType pt = type_of(pos.piece_on(from));
        Bitboard        bb = pt == PAWN ? Attacks::attacks_bb(PAWN, from, us)
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
        return bb;
    }

#if defined(USE_AVX2)
    // exp(-d / 2^FRAC) with 30 fraction bits for 4 values (d <= 32 << FRAC) in 64-bit lanes; identical to
    // expNeg() per lane (d = 32 << FRAC gives 0, like expNeg()).
    static inline __m256i exp_neg4(__m256i d) {
        const __m256i y  = _mm256_srli_epi64(_mm256_mul_epu32(d, _mm256_set1_epi64x(1512775)), FRAC);
        const __m256i nn = _mm256_srli_epi64(y, FRAC);
        const __m256i fr = _mm256_and_si256(y, _mm256_set1_epi64x((i64(1) << FRAC) - 1));
        const __m256i ix = _mm256_srli_epi64(fr, 12);
        const __m256i rm = _mm256_and_si256(fr, _mm256_set1_epi64x(4095));
        // one 64-bit gather of the pair (EXP2_TAB[ix], EXP2_TAB[ix + 1]) instead of two 32-bit gathers
        const __m256i pr = _mm256_i64gather_epi64(reinterpret_cast<const long long*>(g_expPair.v), ix, 8);
        const __m256i t0 = _mm256_and_si256(pr, _mm256_set1_epi64x(0xFFFFFFFF));
        const __m256i t1 = _mm256_srli_epi64(pr, 32);
        __m256i       v  = _mm256_add_epi64(_mm256_mul_epu32(t0, _mm256_sub_epi64(_mm256_set1_epi64x(4096), rm)),
                                            _mm256_mul_epu32(t1, rm));
        return _mm256_srlv_epi64(_mm256_srli_epi64(v, 12), nn);
    }
    static inline i64 hsum_epi64(__m256i v) {
        const __m128i s = _mm_add_epi64(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
        return _mm_cvtsi128_si64(s) + _mm_extract_epi64(s, 1);
    }
    static inline i32 hmax_epi32(__m256i v) {
        __m128i m = _mm_max_epi32(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
        m         = _mm_max_epi32(m, _mm_shuffle_epi32(m, _MM_SHUFFLE(1, 0, 3, 2)));
        m         = _mm_max_epi32(m, _mm_shuffle_epi32(m, _MM_SHUFFLE(2, 3, 0, 1)));
        return _mm_cvtsi128_si32(m);
    }

    // Attention weights of the window of oriented square 'of' (8 entries per vector): for every field
    // logit a = (beta * zhat) / 65536, w = expNeg(max a - a), alpha = (w * (2^60 / sum w)) >> 30.
    // alphaV[slot][k][u] holds alpha of window entry u; eidx[slot][square] maps a square to its entry.
    int window(int of) {
        const int    slot = nslots++;
        const Square from = Square(of ^ flip);
        const uint8_t* drow = g_dir.row[flip ? 1 : 0][from];
        alignas(32) i32 aS[NZ][32 + 8];  // logits a = (beta * zhat) / 65536 (|a| < 2^31), padded
        int             n = 0;
        for (Bitboard b = window_squares(from); b; ++n)
        {
            const Square s = pop_lsb(b);
            const int    d = drow[s];
            eidx[slot][s]  = uint8_t(n);
            aS[0][n]       = i32((g_bq[0][d] * zhv[0][s]) / 65536);
            aS[1][n]       = i32((g_bq[1][d] * zhv[1][s]) / 65536);
        }
        const int     nv   = (n + 7) / 8;
        const __m256i lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
        for (int u = n; u < nv * 8; ++u)
            aS[0][u] = aS[1][u] = 0;  // padding: masked out below

        for (int k = 0; k < NZ; ++k)
        {
            __m256i av[4], valid[4];
            __m256i amaxv = _mm256_set1_epi32(INT32_MIN);
            for (int v = 0; v < nv; ++v)
            {
                av[v]            = _mm256_load_si256(reinterpret_cast<const __m256i*>(&aS[k][8 * v]));
                valid[v]         = _mm256_cmpgt_epi32(_mm256_set1_epi32(n - 8 * v), lane);
                amaxv = _mm256_max_epi32(amaxv, _mm256_blendv_epi8(_mm256_set1_epi32(INT32_MIN), av[v], valid[v]));
            }
            const i32 amax = hmax_epi32(amaxv);

            __m256i wlo[4], whi[4], wsumv = _mm256_setzero_si256();
            for (int v = 0; v < nv; ++v)
            {
                // d = max a - a as unsigned 32 bit, capped at 32 << FRAC (where expNeg() is 0)
                __m256i d = _mm256_sub_epi32(_mm256_set1_epi32(amax), av[v]);
                d         = _mm256_min_epu32(d, _mm256_set1_epi32(32 << FRAC));
                wlo[v] = _mm256_and_si256(exp_neg4(_mm256_cvtepu32_epi64(_mm256_castsi256_si128(d))),
                                          _mm256_cvtepi32_epi64(_mm256_castsi256_si128(valid[v])));
                whi[v] = _mm256_and_si256(exp_neg4(_mm256_cvtepu32_epi64(_mm256_extracti128_si256(d, 1))),
                                          _mm256_cvtepi32_epi64(_mm256_extracti128_si256(valid[v], 1)));
                wsumv  = _mm256_add_epi64(wsumv, _mm256_add_epi64(wlo[v], whi[v]));
            }
            const i64     rinv  = recip_2_60(hsum_epi64(wsumv));
            const __m256i rinvv = _mm256_set1_epi64x(rinv);
            const __m256i pick  = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
            for (int v = 0; v < nv; ++v)
            {
                const __m256i al = _mm256_permutevar8x32_epi32(_mm256_srli_epi64(_mm256_mul_epu32(wlo[v], rinvv), 30), pick);
                const __m256i ah = _mm256_permutevar8x32_epi32(_mm256_srli_epi64(_mm256_mul_epu32(whi[v], rinvv), 30), pick);
                _mm_storeu_si128(reinterpret_cast<__m128i*>(&alphaV[slot][k][8 * v]), _mm256_castsi256_si128(al));
                _mm_storeu_si128(reinterpret_cast<__m128i*>(&alphaV[slot][k][8 * v + 4]), _mm256_castsi256_si128(ah));
            }
        }
        return slot;
    }
#else
    // Portable version of window(); see the AVX2 version above.
    int window(int of) {
        const int    slot = nslots++;
        const Square from = Square(of ^ flip);
        const uint8_t* drow = g_dir.row[flip ? 1 : 0][from];
        i64 a[NZ][32], w[32];
        int n = 0;
        for (Bitboard b = window_squares(from); b; ++n)
        {
            const Square s = pop_lsb(b);
            eidx[slot][s]  = uint8_t(n);
            for (int k = 0; k < NZ; ++k)
                a[k][n] = (g_bq[k][drow[s]] * zhv[k][s]) / 65536;
        }
        for (int k = 0; k < NZ; ++k)
        {
            i64 amax = a[k][0], wsum = 0;
            for (int u = 1; u < n; ++u)
                amax = std::max(amax, a[k][u]);
            for (int u = 0; u < n; ++u)
                wsum += (w[u] = expNeg(amax - a[k][u]));
            const i64 rinv = recip_2_60(wsum);
            for (int u = 0; u < n; ++u)
                alphaV[slot][k][u] = i32((w[u] * rinv) >> 30);
        }
        return slot;
    }
#endif

    // Contribution of the hidden layer to the score of a move of piece type 'type' from oriented square
    // 'of' to the real square 'to' (units 2^-FRAC).
    i64 hidden_score(int of, Square to, int type) {
        int slot = slotMap[of];
        if (slot < 0)
            slotMap[of] = int8_t(slot = window(of));
        const int u = eidx[slot][to];
        return (g_kc[0][type] * alphaV[slot][0][u]) / (i64(1) << 30)
             + (g_kc[1][type] * alphaV[slot][1][u]) / (i64(1) << 30);
    }

    // Field score = bias + sum of the piece-square weights: int8 weights are summed exactly in
    // int32 (at most 32 pieces), then each group's integer multiplier is applied once.
    i64 field(int base, int bias, int gW, int gB, int o, const int16_t* sums) const {
        const int8_t* w   = g_w.data();
        int           acc = 0;
        if (haveSums)
            acc = sums[o];
        else
        {
            const int8_t* row = w + base + o * NF;
            for (int j = 0; j < nfeat; ++j)
                acc += row[fidx[j]];
        }
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
            fF[of] = field(OFF_WF, OFF_BF, G_WF, G_BF, of, sumF), doneF |= u64(1) << of;
        if (!(doneT >> ot & 1))
            fT[ot] = field(OFF_WT, OFF_BT, G_WT, G_BT, ot, sumT), doneT |= u64(1) << ot;

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

        const int eN = theirs.at(to), oN = ours.at(to), eF = theirs.at(from);

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
    AttackCount     ours, theirs;
    int             nfeat = 0;
    u16             fidx[32];
    u64             doneF = 0, doneT = 0;
    i64             fF[64], fT[64];
    bool            haveSums = false;
    // hidden layer state
    Bitboard adj[64] = {};         // piece-mobility graph (real squares)
    alignas(32) i32 zhv[NZ][64];   // zhat (16 fraction bits) per field and real square
    int8_t   slotMap[64];          // oriented from-square -> slot of its attention window (-1: not computed)
    int      nslots = 0;
    alignas(32) i32 alphaV[16][NZ][32];  // attention weights (30 fraction bits): [window slot][field][window entry]
    uint8_t  eidx[16][64];         // [window slot][destination square] -> window entry
#if CRAM_SIMD
    int16_t sumF[64], sumT[64];
#else
    static constexpr const int16_t* sumF = nullptr;
    static constexpr const int16_t* sumT = nullptr;
#endif
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
#if CRAM_SIMD
    std::vector<int8_t>().swap(g_wt);
#endif
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
#if CRAM_SIMD
    g_wt.resize(size_t(2) * 64 * NF);
    for (int f = 0; f < NF; ++f)
        for (int o = 0; o < 64; ++o)
        {
            g_wt[size_t(f) * 64 + o]                   = w[OFF_WF + o * NF + f];
            g_wt[size_t(64) * NF + size_t(f) * 64 + o] = w[OFF_WT + o * NF + f];
        }
#endif
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
    if (g_hidden)
        for (int k = 0; k < NZ; ++k)
        {
            std::fill(g_bq32[k], g_bq32[k] + 16, 0);
            for (int d = 0; d < NDIR; ++d)
            {
                g_bq[k][d]   = i64(g_mult[G_BETA]) * w[OFF_BETA + k * NDIR + d];
                g_bq32[k][d] = i32(g_bq[k][d]);
            }
            for (int t = 0; t < 6; ++t)
                g_kc[k][t] = i64(g_mult[G_KAPPA]) * w[OFF_KAPPA + k * 6 + t];
        }
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

    Ctx ctx(pos, n);
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

    Ctx ctx(pos, n);
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

    Ctx                              ctx(pos, MAX_MOVES);  // not time critical: SIMD path
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
