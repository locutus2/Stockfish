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

// Embed the default net into the binary (same mechanism as the NNUE net, see network.cpp).
// Unlike the NNUE net, this does NOT depend on NNUE_EMBEDDING_OFF: Fishtest workers build with
// CXXFLAGS=-DNNUE_EMBEDDING_OFF (they copy the NNUE net next to the binary), but they do not
// provide the CRAM net, so it must always be compiled in. MSVC falls back to loading the file.
#define INCBIN_SILENCE_BITCODE_WARNING
#include "incbin/incbin.h"

#if !defined(UNIVERSAL_BINARY) && !defined(_MSC_VER)
INCBIN(EmbeddedCRAM, "model-i8.bin");
#elif defined(UNIVERSAL_BINARY)
// Provided by universal/cram_embed.cpp
extern const unsigned char gEmbeddedCRAMData[];
extern const unsigned int  gEmbeddedCRAMSize;
#else
const unsigned char gEmbeddedCRAMData[1] = {0x0};
const unsigned int  gEmbeddedCRAMSize    = 1;
#endif

namespace Stockfish::CRAM {

bool g_loaded = false;
int  g_plies  = 256;

namespace {

std::vector<int8_t> g_w;  // the net, NPARAM int8 weights; only replaced while no search runs
i32                 g_mult[NGROUPS];   // value of weight i of group g = g_w[i] * g_mult[g] / 2^FRAC
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

        return fF[of] + fT[ot] + i64(g_mult[G_K]) * w[OFF_K + of * 64 + ot]
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

void unload() {
    g_loaded = false;
    std::vector<int8_t>().swap(g_w);
}

void set_weight(int w) { g_weight = w; }
void set_plies(int p) { g_plies = p; }
void set_gate(int pct) { g_gate = pct; }
int  gate_percent() { return g_gate; }

namespace {

// Parses a CRM3 net from memory and installs it. Returns an error message, or an empty string.
std::string parse_net(const unsigned char* data, size_t size) {
    constexpr size_t HEADER = 4 * sizeof(u32) + NGROUPS * sizeof(i32);

    if (size < 2 * sizeof(u32))
        return "cannot read header";

    u32 hdr[4] = {0, 0, 0, 0};
    std::memcpy(hdr, data, 2 * sizeof(u32));  // little-endian host assumed, like the NNUE loader
    if (hdr[0] == MAGIC_F32 || hdr[0] == MAGIC_Q8F)
        return "old net with float weights/scales (CRM1/CRM2) is no longer supported; "
               "retrain with the integer chess_cram";
    if (size < HEADER)
        return "wrong file size";
    std::memcpy(hdr, data, 4 * sizeof(u32));
    if (hdr[0] != MAGIC || hdr[1] != u32(NPARAM) || hdr[2] != u32(NGROUPS) || hdr[3] != u32(FRAC))
        return "not a compatible int8 net (bad magic, parameter count, group count or "
               "fixed-point format; expected "
             + std::to_string(NPARAM) + " parameters)";
    if (size != HEADER + size_t(NPARAM))
        return "wrong file size";

    i32 mult[NGROUPS];
    std::memcpy(mult, data + 4 * sizeof(u32), sizeof(mult));
    for (i32 v : mult)
        if (v < 1 || v > MAX_MULT)
            return "invalid multiplier";

    const int8_t* w = reinterpret_cast<const int8_t*>(data + HEADER);
    g_w.assign(w, w + NPARAM);
    std::copy(mult, mult + NGROUPS, g_mult);
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
        return "CRAM int8 net loaded: embedded " + std::string(DefaultFile);
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
    return "CRAM int8 net loaded: " + path;
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
