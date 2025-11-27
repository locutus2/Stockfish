/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2025 The Stockfish developers (see AUTHORS file)

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "movepick.h"

#include <cassert>
#include <limits>
#include <utility>

#include "bitboard.h"
#include "misc.h"
#include "position.h"

namespace Stockfish {

namespace {

enum Stages {
    // generate main search moves
    MAIN_TT,
    CAPTURE_INIT,
    GOOD_CAPTURE,
    QUIET_INIT,
    GOOD_QUIET,
    BAD_CAPTURE,
    BAD_QUIET,

    // generate evasion moves
    EVASION_TT,
    EVASION_INIT,
    EVASION,

    // generate probcut moves
    PROBCUT_TT,
    PROBCUT_INIT,
    PROBCUT,

    // generate qsearch moves
    QSEARCH_TT,
    QCAPTURE_INIT,
    QCAPTURE
};


// Sort moves in descending order up to and including a given limit.
// The order of moves smaller than the limit is left unspecified.
void partial_insertion_sort(ExtMove* begin, ExtMove* end, int limit) {

    for (ExtMove *sortedEnd = begin, *p = begin + 1; p < end; ++p)
        if (p->value >= limit)
        {
            ExtMove tmp = *p, *q;
            *p          = *++sortedEnd;
            for (q = sortedEnd; q != begin && *(q - 1) < tmp; --q)
                *q = *(q - 1);
            *q = tmp;
        }
}

}  // namespace


// Constructors of the MovePicker class. As arguments, we pass information
// to decide which class of moves to emit, to help sorting the (presumably)
// good moves first, and how important move ordering is at the current node.

// MovePicker constructor for the main search and for the quiescence search
MovePicker::MovePicker(const Position&              p,
                       Move                         ttm,
                       Depth                        d,
                       const ButterflyHistory*      mh,
                       const LowPlyHistory*         lph,
                       const CapturePieceToHistory* cph,
                       const PieceToHistory**       ch,
                       const PawnHistory*           ph,
                       int                          pl) :
    pos(p),
    mainHistory(mh),
    lowPlyHistory(lph),
    captureHistory(cph),
    continuationHistory(ch),
    pawnHistory(ph),
    ttMove(ttm),
    depth(d),
    ply(pl) {

    if (pos.checkers())
        stage = EVASION_TT + !(ttm && pos.pseudo_legal(ttm));

    else
        stage = (depth > 0 ? MAIN_TT : QSEARCH_TT) + !(ttm && pos.pseudo_legal(ttm));
}

// MovePicker constructor for ProbCut: we generate captures with Static Exchange
// Evaluation (SEE) greater than or equal to the given threshold.
MovePicker::MovePicker(const Position& p, Move ttm, int th, const CapturePieceToHistory* cph) :
    pos(p),
    captureHistory(cph),
    ttMove(ttm),
    threshold(th) {
    assert(!pos.checkers());

    stage = PROBCUT_TT + !(ttm && pos.capture_stage(ttm) && pos.pseudo_legal(ttm));
}

// Assigns a numerical value to each move in a list, used for sorting.
// Captures are ordered by Most Valuable Victim (MVV), preferring captures
// with a good history. Quiets moves are ordered using the history tables.
template<GenType Type>
ExtMove* MovePicker::score(MoveList<Type>& ml) {

    static_assert(Type == CAPTURES || Type == QUIETS || Type == EVASIONS, "Wrong type");

    Color us = pos.side_to_move();

    [[maybe_unused]] Bitboard threatByLesser[KING + 1];
    if constexpr (Type == QUIETS)
    {
        threatByLesser[PAWN]   = 0;
        threatByLesser[KNIGHT] = threatByLesser[BISHOP] = pos.attacks_by<PAWN>(~us);
        threatByLesser[ROOK] =
          pos.attacks_by<KNIGHT>(~us) | pos.attacks_by<BISHOP>(~us) | threatByLesser[KNIGHT];
        threatByLesser[QUEEN] = pos.attacks_by<ROOK>(~us) | threatByLesser[ROOK];
        threatByLesser[KING]  = pos.attacks_by<QUEEN>(~us) | threatByLesser[QUEEN];
    }

    ExtMove* it = cur;
    for (auto move : ml)
    {
        ExtMove& m = *it++;
        m          = move;

        const Square    from          = m.from_sq();
        const Square    to            = m.to_sq();
        const Piece     pc            = pos.moved_piece(m);
        const PieceType pt            = type_of(pc);
        const Piece     capturedPiece = pos.piece_on(to);

        if constexpr (Type == CAPTURES)
            m.value = (*captureHistory)[pc][to][type_of(capturedPiece)]
                    + 7 * int(PieceValue[capturedPiece]);

        else if constexpr (Type == QUIETS)
        {
            // histories
            m.value = 2 * (*mainHistory)[us][m.raw()];
            m.value += 2 * (*pawnHistory)[pawn_history_index(pos)][pc][to];
            m.value += (3 * (*continuationHistory[0])[pc][to] + (*continuationHistory[2])[pc][to]) / 2;
            m.value += (*continuationHistory[1])[pc][to];
            m.value += (*continuationHistory[3])[pc][to];
            m.value += (*continuationHistory[5])[pc][to];
            std::vector<int> x = {  
                (*mainHistory)[us][m.raw()],
                (*pawnHistory)[pawn_history_index(pos)][pc][to],
                (*continuationHistory[0])[pc][to],
                (*continuationHistory[1])[pc][to],
                (*continuationHistory[2])[pc][to],
                (*continuationHistory[3])[pc][to],
                (*continuationHistory[5])[pc][to],
            };

              std::vector<int> h = {m.value};
              //offset: 4.16272 factors: 0.000333064,0.000277362,0.000118239,0.000168658,0.000124153,0.00018503,0.000185863
              double offset = 4.16272;
              constexpr int SCALE = 1000;
              std::vector<double> factor = { 0.000333064,0.000277362,0.000118239,0.000168658,0.000124153,0.00018503,0.000185863};
              h.push_back(((factor[0] * x[0]
                  + factor[1] * x[1]
                  + factor[2] * x[2]
                  + factor[3] * x[3]
                  + factor[4] * x[4]
                  + factor[5] * x[5]
                  + factor[6] * x[6]
                  + offset) * SCALE - 0.039791) * 6.4600191420864874271295571217263 - 26316.2);
              /*
               * Mean #1000: Total 1009242979 Mean -26316.2
               * Mean #1001: Total 1009242979 Mean 0.039791
               * Stdev #1000: Total 1009242979 Stdev 20788.6
               * Stdev #1001: Total 1009242979 Stdev 3218.04
               * Correl. #1000: Total 1009242979 Coefficient 0.992495 Cov 6.63966e+07
               * */
              /*
              h.push_back((0.553 * x[0]
                  + 0.465 * x[1]
                  + 0.847 * x[2]
                  + 0.991 * x[3]
                  + 0.881 * x[4]
                  + 1.043 * x[5]
                  + 1 * x[6]
                  +15333.2) / 16062.2 * 20788.6 / 22334.3 * 20788.6 - 26316.2 + 884.5);
                  */

              dbg_mean_of(h[0], 1000);
              dbg_stdev_of(h[0], 1000);
              dbg_mean_of(h[1], 1001);
              dbg_stdev_of(h[1], 1001);
              dbg_correl_of(h[0], h[1], 1000);

              //PCA::add(x);
              //m.value = h[1];
              /*
               * Hit #0: Total 116532149 Hits 13600045 Hit Rate (%) 11.6706
               * Mean #0: Total 116532149 Mean -1745.23
               * Mean #1: Total 116532149 Mean -285.684
               * Mean #2: Total 116532149 Mean 1524.99
               * Mean #3: Total 116532149 Mean -528.671
               * Mean #4: Total 116532149 Mean -640.566
               * Mean #5: Total 116532149 Mean -1384.11
               * Mean #6: Total 116532149 Mean -1198.41
               * Mean #7: Total 116532149 Mean -505.138
               * Mean #100: Total 1009242979 Mean -26316.2
               * Mean #101: Total 1009242979 Mean -15333.2
               * Stdev #0: Total 116532149 Stdev 4038.44
               * Stdev #1: Total 116532149 Stdev 4203.3
               * Stdev #2: Total 116532149 Stdev 8193
               * Stdev #3: Total 116532149 Stdev 6698.41
               * Stdev #4: Total 116532149 Stdev 7619.97
               * Stdev #5: Total 116532149 Stdev 6675.44
               * Stdev #6: Total 116532149 Stdev 6384.15
               * Stdev #7: Total 116532149 Stdev 3127.28
               * Stdev #100: Total 1009242979 Stdev 20788.6
               * Stdev #101: Total 1009242979 Stdev 16062.2
               * Correl. #0: Total 116532149 Coefficient -0.2085 Cov -0.426611
               * Correl. #10: Total 116532149 Coefficient 0.212801 Cov 275.923
               * Correl. #11: Total 116532149 Coefficient 0.194474 Cov 262.453
               * Correl. #12: Total 116532149 Coefficient 0.198944 Cov 523.327
               * Correl. #13: Total 116532149 Coefficient 0.141724 Cov 304.801
               * Correl. #14: Total 116532149 Coefficient 0.122994 Cov 300.909
               * Correl. #15: Total 116532149 Coefficient 0.162257 Cov 347.763
               * Correl. #16: Total 116532149 Coefficient 0.139336 Cov 285.606
               * Correl. #17: Total 116532149 Coefficient 0.275038 Cov 276.16
               * Correl. #100: Total 1009242979 Coefficient 0.931169 Cov 3.10928e+08
               *
               * Hit #0: Total 116532149 Hits 13600045 Hit Rate (%) 11.6706
               * Mean #0: Total 116532149 Mean -1745.23
               * Mean #1: Total 116532149 Mean -285.684
               * Mean #2: Total 116532149 Mean 1524.99
               * Mean #3: Total 116532149 Mean -528.671
               * Mean #4: Total 116532149 Mean -640.566
               * Mean #5: Total 116532149 Mean -1384.11
               * Mean #6: Total 116532149 Mean -1198.41
               * Mean #7: Total 116532149 Mean -505.138
               * Mean #100: Total 1009242979 Mean -26316.2
               * Mean #101: Total 1009242979 Mean -26316.3
               * Stdev #0: Total 116532149 Stdev 4038.44
               * Stdev #1: Total 116532149 Stdev 4203.3
               * Stdev #2: Total 116532149 Stdev 8193
               * Stdev #3: Total 116532149 Stdev 6698.41
               * Stdev #4: Total 116532149 Stdev 7619.97
               * Stdev #5: Total 116532149 Stdev 6675.44
               * Stdev #6: Total 116532149 Stdev 6384.15
               * Stdev #7: Total 116532149 Stdev 3127.28
               * Stdev #100: Total 1009242979 Stdev 20788.6
               * Stdev #101: Total 1009242979 Stdev 20788.8
               * Correl. #0: Total 116532149 Coefficient -0.2085 Cov -0.426611
               * Correl. #10: Total 116532149 Coefficient 0.212801 Cov 275.923
               * Correl. #11: Total 116532149 Coefficient 0.194474 Cov 262.453
               * Correl. #12: Total 116532149 Coefficient 0.198944 Cov 523.327
               * Correl. #13: Total 116532149 Coefficient 0.141724 Cov 304.801
               * Correl. #14: Total 116532149 Coefficient 0.122994 Cov 300.909
               * Correl. #15: Total 116532149 Coefficient 0.162257 Cov 347.763
               * Correl. #16: Total 116532149 Coefficient 0.139336 Cov 285.606
               * Correl. #17: Total 116532149 Coefficient 0.275038 Cov 276.16
               * Correl. #100: Total 1009242979 Coefficient 0.931169 Cov 4.02423e+08
               * PCA N=0
               *
               * ===========================
               * Total time (ms) : 645724
               * Nodes searched  : 191897230
               * Nodes/second    : 297181
               *
               //m.value = h[1];
               * Hit #0: Total 126900978 Hits 14762863 Hit Rate (%) 11.6334
                                   * Mean #0: Total 126900978 Mean -1674.42
                                   * Mean #1: Total 126900978 Mean -372.984
                                   * Mean #2: Total 126900978 Mean 1269.26
                                   * Mean #3: Total 126900978 Mean -684.834
                                   * Mean #4: Total 126900978 Mean -632.925
                                   * Mean #5: Total 126900978 Mean -1367.07
                                   * Mean #6: Total 126900978 Mean -1208.68
                                   * Mean #7: Total 126900978 Mean -556.377
                                   * Mean #100: Total 1098838922 Mean -26673.7
                                   * Mean #101: Total 1098838922 Mean -27200.7
                                   * Stdev #0: Total 126900978 Stdev 4107.02
                                   * Stdev #1: Total 126900978 Stdev 4282.69
                                   * Stdev #2: Total 126900978 Stdev 8587.97
                                   * Stdev #3: Total 126900978 Stdev 7080.77
                                   * Stdev #4: Total 126900978 Stdev 8014.07
                                   * Stdev #5: Total 126900978 Stdev 7062.77
                                   * Stdev #6: Total 126900978 Stdev 6772.59
                                   * Stdev #7: Total 126900978 Stdev 3622.03
                                   * Stdev #100: Total 1098838922 Stdev 22606.5
                                   * Stdev #101: Total 1098838922 Stdev 22334.3
                                   * Correl. #0: Total 126900978 Coefficient -0.209791 Cov -0.436185
                                   * Correl. #10: Total 126900978 Coefficient 0.211742 Cov 278.825
                                   * Correl. #11: Total 126900978 Coefficient 0.196021 Cov 269.164
                                   * Correl. #12: Total 126900978 Coefficient 0.218513 Cov 601.679
                                   * Correl. #13: Total 126900978 Coefficient 0.172218 Cov 390.981
                                   * Correl. #14: Total 126900978 Coefficient 0.147643 Cov 379.371
                                   * Correl. #15: Total 126900978 Coefficient 0.186845 Cov 423.11
                                   * Correl. #16: Total 126900978 Coefficient 0.166587 Cov 361.736
                                   * Correl. #17: Total 126900978 Coefficient 0.284583 Cov 330.491
                                   * Correl. #100: Total 1098838922 Coefficient 0.940368 Cov 4.74792e+08
                                   * PCA N=0
                                   *
                                   * ===========================
                                   * Total time (ms) : 716244
                                                       * Nodes searched  : 209467861
                                                       * Nodes/second    : 292453
                                                       *
               * Hit #0: Total 123020115 Hits 14352652 Hit Rate (%) 11.6669
                                                       * Mean #0: Total 123020115 Mean -1721.84
                                                       * Mean #1: Total 123020115 Mean -387.424
                                                       * Mean #2: Total 123020115 Mean 1370.85
                                                       * Mean #3: Total 123020115 Mean -612.142
                                                       * Mean #4: Total 123020115 Mean -552.831
                                                       * Mean #5: Total 123020115 Mean -1289.6
                                                       * Mean #6: Total 123020115 Mean -1147.05
                                                       * Mean #7: Total 123020115 Mean -508.073
                                                       * Mean #100: Total 1064161878 Mean -26397.9
                                                       * Mean #101: Total 1064161878 Mean -25909.8
                                                       * Stdev #0: Total 123020115 Stdev 4086.1
                                                       * Stdev #1: Total 123020115 Stdev 4265.47
                                                       * Stdev #2: Total 123020115 Stdev 8312.71
                                                       * Stdev #3: Total 123020115 Stdev 6816.67
                                                       * Stdev #4: Total 123020115 Stdev 7728.61
                                                       * Stdev #5: Total 123020115 Stdev 6793.75
                                                       * Stdev #6: Total 123020115 Stdev 6508.08
                                                       * Stdev #7: Total 123020115 Stdev 3384.48
                                                       * Stdev #100: Total 1064161878 Stdev 21748.8
                                                       * Stdev #101: Total 1064161878 Stdev 19735
                                                       * Correl. #0: Total 123020115 Coefficient -0.210344 Cov -0.434365
                                                       * Correl. #10: Total 123020115 Coefficient 0.209924 Cov 275.367
                                                       * Correl. #11: Total 123020115 Coefficient 0.193397 Cov 264.823
                                                       * Correl. #12: Total 123020115 Coefficient 0.214296 Cov 571.869
                                                       * Correl. #13: Total 123020115 Coefficient 0.162407 Cov 355.4
                                                       * Correl. #14: Total 123020115 Coefficient 0.137889 Cov 342.113
                                                       * Correl. #15: Total 123020115 Coefficient 0.177521 Cov 387.168
                                                       * Correl. #16: Total 123020115 Coefficient 0.156025 Cov 325.976
                                                       * Correl. #17: Total 123020115 Coefficient 0.281758 Cov 306.132
                                                       * Correl. #100: Total 1064161878 Coefficient 0.936585 Cov 4.01994e+08
                                                       * PCA N=0
                                                       *
                                                       * ===========================
                                                       * Total time (ms) : 723877
                                                                           * Nodes searched  : 203259689
                                                                           * Nodes/second    : 280793
               */

            // bonus for checks
            m.value += (bool(pos.check_squares(pt) & to) && pos.see_ge(m, -75)) * 16384;

            // penalty for moving to a square threatened by a lesser piece
            // or bonus for escaping an attack by a lesser piece.
            int v = threatByLesser[pt] & to ? -19 : 20 * bool(threatByLesser[pt] & from);
            m.value += PieceValue[pt] * v;


            if (ply < LOW_PLY_HISTORY_SIZE)
                m.value += 8 * (*lowPlyHistory)[ply][m.raw()] / (1 + ply);
        }

        else  // Type == EVASIONS
        {
            if (pos.capture_stage(m))
                m.value = PieceValue[capturedPiece] + (1 << 28);
            else
            {
                m.value = (*mainHistory)[us][m.raw()] + (*continuationHistory[0])[pc][to];
                if (ply < LOW_PLY_HISTORY_SIZE)
                    m.value += (*lowPlyHistory)[ply][m.raw()];
            }
        }
    }
    return it;
}

// Returns the next move satisfying a predicate function.
// This never returns the TT move, as it was emitted before.
template<typename Pred>
Move MovePicker::select(Pred filter) {

    for (; cur < endCur; ++cur)
        if (*cur != ttMove && filter())
            return *cur++;

    return Move::none();
}

// This is the most important method of the MovePicker class. We emit one
// new pseudo-legal move on every call until there are no more moves left,
// picking the move with the highest score from a list of generated moves.
Move MovePicker::next_move() {

    constexpr int goodQuietThreshold = -14000;
top:
    switch (stage)
    {

    case MAIN_TT :
    case EVASION_TT :
    case QSEARCH_TT :
    case PROBCUT_TT :
        ++stage;
        return ttMove;

    case CAPTURE_INIT :
    case PROBCUT_INIT :
    case QCAPTURE_INIT : {
        MoveList<CAPTURES> ml(pos);

        cur = endBadCaptures = moves;
        endCur = endCaptures = score<CAPTURES>(ml);

        partial_insertion_sort(cur, endCur, std::numeric_limits<int>::min());
        ++stage;
        goto top;
    }

    case GOOD_CAPTURE :
        if (select([&]() {
                if (pos.see_ge(*cur, -cur->value / 18))
                    return true;
                std::swap(*endBadCaptures++, *cur);
                return false;
            }))
            return *(cur - 1);

        ++stage;
        [[fallthrough]];

    case QUIET_INIT :
        if (!skipQuiets)
        {
            MoveList<QUIETS> ml(pos);

            endCur = endGenerated = score<QUIETS>(ml);

            partial_insertion_sort(cur, endCur, -3560 * depth);
        }

        ++stage;
        [[fallthrough]];

    case GOOD_QUIET :
        if (!skipQuiets && select([&]() { return cur->value > goodQuietThreshold; }))
            return *(cur - 1);

        // Prepare the pointers to loop over the bad captures
        cur    = moves;
        endCur = endBadCaptures;

        ++stage;
        [[fallthrough]];

    case BAD_CAPTURE :
        if (select([]() { return true; }))
            return *(cur - 1);

        // Prepare the pointers to loop over quiets again
        cur    = endCaptures;
        endCur = endGenerated;

        ++stage;
        [[fallthrough]];

    case BAD_QUIET :
        if (!skipQuiets)
            return select([&]() { return cur->value <= goodQuietThreshold; });

        return Move::none();

    case EVASION_INIT : {
        MoveList<EVASIONS> ml(pos);

        cur    = moves;
        endCur = endGenerated = score<EVASIONS>(ml);

        partial_insertion_sort(cur, endCur, std::numeric_limits<int>::min());
        ++stage;
        [[fallthrough]];
    }

    case EVASION :
    case QCAPTURE :
        return select([]() { return true; });

    case PROBCUT :
        return select([&]() { return pos.see_ge(*cur, threshold); });
    }

    assert(false);
    return Move::none();  // Silence warning
}

void MovePicker::skip_quiet_moves() { skipQuiets = true; }

}  // namespace Stockfish
