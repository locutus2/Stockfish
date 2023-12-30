/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2023 The Stockfish developers (see AUTHORS file)

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

#include <algorithm>
#include <cassert>
#include <iterator>
#include <utility>
#include <vector>

#include "bitboard.h"
#include "position.h"
#include "stats.h"
#include "thread.h"

namespace Stockfish {

int HISTORY_WEIGHT[N_HISTORY][2];

int Dmax, Dmin;
//constexpr int EvCapOffset = 45334;
//constexpr int EvCapWeight = 1;

//constexpr int EvCapOffset = -16778;
//constexpr int EvCapWeight = 12;

constexpr int EvCapOffset = QueenValue;
//constexpr int EvCapOffset = (1 << 28);
constexpr int EvCapWeight = 1;

void calculateMinMax(StageType stage, int& Dmin, int& Dmax);

void calculateMinMax(StageType stage, int& Min, int& Max) {
    for (int i = 0; i < N_HISTORY; ++i)
    {
        int S = HISTORY_WEIGHT_MASTER[stage][i][1] * HISTORY_WEIGHT[i][1];
        int W = HISTORY_WEIGHT_MASTER[stage][i][0] * HISTORY_WEIGHT[i][1]
              + HISTORY_WEIGHT[i][0] * HISTORY_WEIGHT_MASTER[stage][i][1];
        Max += HISTORY_RANGE[i][W > 0] * W / S;
        Min += HISTORY_RANGE[i][W <= 0] * W / S;
    }
}

void init_stats(bool onlyD) {

    if (!onlyD)
    {
        for (int i = 0; i < N_HISTORY; ++i)
        {
            HISTORY_WEIGHT[i][0] = HISTORY_START[i][0];
            HISTORY_WEIGHT[i][1] = HISTORY_START[i][1];
        }
    }

    Dmax = 0;
    Dmin = 0;

    if (STATS_QUIETS)
    {
        calculateMinMax(STAGE_QUIETS, Dmin, Dmax);
    }
    else if (STATS_REFUTATION)
    {
        Dmax = 1;
        Dmax = -1;
        calculateMinMax(STAGE_REFUTATIONS, Dmin, Dmax);
    }
    else if (STATS_CAPTURE_MAIN)
    {
        calculateMinMax(STAGE_CAPTURES_MAIN, Dmin, Dmax);
        Dmax += 7 * QueenValue / 16;
    }
    else if (STATS_QUIET_EVASION_MAIN && STATS_CAPTURE_EVASION_MAIN)
    {
        int DmaxCap = 0;
        int DminCap = 0;
        calculateMinMax(STAGE_CAPTURE_EVASIONS, DminCap, DmaxCap);
        DmaxCap += EvCapWeight * QueenValue;
        DminCap -= EvCapWeight * int(KING);
        DmaxCap += EvCapOffset;
        DminCap += EvCapOffset;

        int DmaxQuiet = 0;
        int DminQuiet = 0;
        calculateMinMax(STAGE_QUIET_EVASIONS, DminQuiet, DmaxQuiet);

        Dmax = std::max(DmaxCap, DmaxQuiet);
        Dmin = std::min(DminCap, DminQuiet);
    }
    else if (STATS_QUIET_EVASION_MAIN || STATS_QUIET_EVASION_QS)
    {
        calculateMinMax(STAGE_QUIET_EVASIONS, Dmin, Dmax);
    }
    else if (STATS_CAPTURE_EVASION_MAIN || STATS_CAPTURE_EVASION_QS)
    {
        calculateMinMax(STAGE_CAPTURE_EVASIONS, Dmin, Dmax);
        Dmax += EvCapWeight * QueenValue;
        Dmin -= EvCapWeight * int(KING);
        Dmax += EvCapOffset;
        Dmin += EvCapOffset;
    }

    Dmax = std::max(Dmax, Dmin + 1);
    //std::cerr << "[" << Dmin << "|" << Dmax << "]" << std::flush << std::endl;
    //std::exit(1);
}

namespace {

enum Stages {
    // generate main search moves
    MAIN_TT,
    CAPTURE_INIT,
    GOOD_CAPTURE,
    REFUTATION,
    QUIET_INIT,
    QUIET,
    BAD_CAPTURE,

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
    QCAPTURE,
    QCHECK_INIT,
    QCHECK
};

// Sort moves in descending order up to and including
// a given limit. The order of moves smaller than the limit is left unspecified.
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
// to help it return the (presumably) good moves first, to decide which
// moves to return (in the quiescence search, for instance, we only want to
// search captures, promotions, and some checks) and how important a good
// move ordering is at the current node.

// MovePicker constructor for the main search
MovePicker::MovePicker(const Position&              p,
                       Move                         ttm,
                       Depth                        d,
                       const ButterflyHistory*      mh,
                       const CapturePieceToHistory* cph,
                       const PieceToHistory**       ch,
                       const PawnHistory*           ph,
                       //const PieceToHistory***       kh,
                       Move        cm,
                       const Move* killers,
                       bool        cond) :
    pos(p),
    mainHistory(mh),
    captureHistory(cph),
    continuationHistory(ch),
    pawnHistory(ph),
    //kingHistory(kh),
    ttMove(ttm),
    refutations{
      {killers[0], 1 },
      {killers[1], 0 },
      {cm,         -1}
},
    depth(d), C(cond) {
    assert(d > 0);

    stage = (pos.checkers() ? EVASION_TT : MAIN_TT) + !(ttm && pos.pseudo_legal(ttm));
}

// Constructor for quiescence search
MovePicker::MovePicker(const Position&              p,
                       Move                         ttm,
                       Depth                        d,
                       const ButterflyHistory*      mh,
                       const CapturePieceToHistory* cph,
                       const PieceToHistory**       ch,
                       const PawnHistory*           ph,
                       bool                         cond) :
    pos(p),
    mainHistory(mh),
    captureHistory(cph),
    continuationHistory(ch),
    pawnHistory(ph),
    ttMove(ttm),
    depth(d),
    C(cond) {
    assert(d <= 0);

    stage = (pos.checkers() ? EVASION_TT : QSEARCH_TT) + !(ttm && pos.pseudo_legal(ttm));
}

// Constructor for ProbCut: we generate captures with SEE greater
// than or equal to the given threshold.
MovePicker::MovePicker(const Position& p, Move ttm, Value th, const CapturePieceToHistory* cph) :
    pos(p),
    captureHistory(cph),
    ttMove(ttm),
    threshold(th),
    C(false) {
    assert(!pos.checkers());

    stage = PROBCUT_TT
          + !(ttm && pos.capture_stage(ttm) && pos.pseudo_legal(ttm) && pos.see_ge(ttm, threshold));
}

// Assigns a numerical value to each move in a list, used
// for sorting. Captures are ordered by Most Valuable Victim (MVV), preferring
// captures with a good history. Quiets moves are ordered using the history tables.
template<ScoreType Type>
void MovePicker::score() {

    static_assert(Type == SCORE_CAPTURES || Type == SCORE_REFUTATIONS || Type == SCORE_QUIETS
                    || Type == SCORE_EVASIONS,
                  "Wrong type");

    [[maybe_unused]] Bitboard threatenedByPawn, threatenedByMinor, threatenedByRook,
      threatenedPieces;
    if constexpr (Type == SCORE_QUIETS)
    {
        Color us = pos.side_to_move();

        threatenedByPawn = pos.attacks_by<PAWN>(~us);
        threatenedByMinor =
          pos.attacks_by<KNIGHT>(~us) | pos.attacks_by<BISHOP>(~us) | threatenedByPawn;
        threatenedByRook = pos.attacks_by<ROOK>(~us) | threatenedByMinor;

        // Pieces threatened by pieces of lesser material value
        threatenedPieces = (pos.pieces(us, QUEEN) & threatenedByRook)
                         | (pos.pieces(us, ROOK) & threatenedByMinor)
                         | (pos.pieces(us, KNIGHT, BISHOP) & threatenedByPawn);
    }

    int           k          = 0;
    constexpr int position[] = {4096, 0, -4096};
    //constexpr int position[] = {1 - 7183 / 8, -7183 / 8, 7183 / 8};

    for (auto& m : *this)
    {
        if constexpr (Type == SCORE_CAPTURES)
            m.value =
              (7 * int(PieceValue[pos.piece_on(to_sq(m))])
               + (*captureHistory)[pos.moved_piece(m)][to_sq(m)][type_of(pos.piece_on(to_sq(m)))])
              / 16;

        else if constexpr (Type == SCORE_REFUTATIONS)
        {}

        else if constexpr (Type == SCORE_QUIETS)
        {
            Piece     pc   = pos.moved_piece(m);
            PieceType pt   = type_of(pos.moved_piece(m));
            Square    from = from_sq(m);
            Square    to   = to_sq(m);

            // histories
            m.value = 2 * (*mainHistory)[pos.side_to_move()][from_to(m)];
            m.value += 2 * (*pawnHistory)[pawn_structure(pos)][pc][to];
            m.value += 2 * (*continuationHistory[0])[pc][to];
            m.value += (*continuationHistory[1])[pc][to];
            m.value += (*continuationHistory[2])[pc][to] / 4;
            m.value += (*continuationHistory[3])[pc][to];
            m.value += (*continuationHistory[5])[pc][to];
            //int v = m.value;
            // bonus for checks
            m.value += bool(pos.check_squares(pt) & to) * 16384;

            // bonus for escaping from capture
            m.value += threatenedPieces & from ? (pt == QUEEN && !(to & threatenedByRook)   ? 50000
                                                  : pt == ROOK && !(to & threatenedByMinor) ? 25000
                                                  : !(to & threatenedByPawn)                ? 15000
                                                                                            : 0)
                                               : 0;

            // malus for putting piece en prise
            m.value -= !(threatenedPieces & from)
                       ? (pt == QUEEN ? bool(to & threatenedByRook) * 50000
                                          + bool(to & threatenedByMinor) * 10000
                                          + bool(to & threatenedByPawn) * 20000
                          : pt == ROOK ? bool(to & threatenedByMinor) * 25000
                                           + bool(to & threatenedByPawn) * 10000
                          : pt != PAWN ? bool(to & threatenedByPawn) * 15000
                                       : 0)
                       : 0;
        }

        else  // Type == SCORE_EVASIONS
        {
            if (pos.capture_stage(m))
            {
                m.value =
                  EvCapWeight
                    * (PieceValue[pos.piece_on(to_sq(m))] - Value(type_of(pos.moved_piece(m))))
                  + EvCapOffset;
                //dbg_mean_of(m.value - EvCapOffset, 100);
                //dbg_stdev_of(m.value - EvCapOffset, 100);
            }
            else
            {
                m.value = (*mainHistory)[pos.side_to_move()][from_to(m)]
                        + (*continuationHistory[0])[pos.moved_piece(m)][to_sq(m)]
                        + (*pawnHistory)[pawn_structure(pos)][pos.moved_piece(m)][to_sq(m)];
                //dbg_mean_of(m.value, 101);
                //dbg_stdev_of(m.value, 101);
            }
        }

        if (C
            && ((STATS_QUIETS && Type == SCORE_QUIETS)
                || (STATS_REFUTATION && Type == SCORE_REFUTATIONS)
                || (STATS_QUIET_EVASION_MAIN && Type == SCORE_EVASIONS && depth > 0)
                || (STATS_QUIET_EVASION_QS && Type == SCORE_EVASIONS && depth <= 0)
                || (STATS_CAPTURE_EVASION_MAIN && Type == SCORE_EVASIONS && depth > 0)
                || (STATS_CAPTURE_EVASION_QS && Type == SCORE_EVASIONS && depth <= 0)
                || (STATS_CAPTURE_MAIN && Type == SCORE_CAPTURES && depth > 0)))
        {
            Piece     pc   = pos.moved_piece(m);
            PieceType pt   = type_of(pos.moved_piece(m));
            Square    from = from_sq(m);
            Square    to   = to_sq(m);

            int V                 = 0;
            int values[N_HISTORY] = {
              (*mainHistory)[pos.side_to_move()][from_to(m)],
              (*pawnHistory)[pawn_structure(pos)][pc][to],
              pos.this_thread()->inCheckHistory[pos.side_to_move()][from_to(m)],
              (*continuationHistory[0])[pc][to],
              (*continuationHistory[1])[pc][to],
              (*continuationHistory[2])[pc][to],
              (*continuationHistory[3])[pc][to],
              (*continuationHistory[4])[pc][to],
              (*continuationHistory[5])[pc][to],
              (*captureHistory)[pc][to][type_of(pos.piece_on(to_sq(m)))],
              std::max(int((*continuationHistory[0])[pc][to]), 0),
              std::min(int((*continuationHistory[0])[pc][to]), 0),
              (*mainHistory)[pos.side_to_move()][from_to(m)]
                * (*pawnHistory)[pawn_structure(pos)][pc][to] / 8192,
              (*mainHistory)[pos.side_to_move()][from_to(m)]
                * (8192 + (*pawnHistory)[pawn_structure(pos)][pc][to]) / (2 * 8192),
              (7183 + (*mainHistory)[pos.side_to_move()][from_to(m)])
                * (8192 + (*pawnHistory)[pawn_structure(pos)][pc][to]) / (2 * 8192),
              (STATS_REFUTATION ? position[k] : 0),
              bool(pos.check_squares(pt) & to) * 16384,
              (STATS_QUIETS && threatenedPieces & from
                 ? (pt == QUEEN && !(to & threatenedByRook)   ? 50000
                    : pt == ROOK && !(to & threatenedByMinor) ? 25000
                    : !(to & threatenedByPawn)                ? 15000
                                                              : 0)
                 : 0),
              (STATS_QUIETS && !(threatenedPieces & from)
                 ? (pt == QUEEN
                      ? bool(to & threatenedByRook) * 50000 + bool(to & threatenedByMinor) * 10000
                          + bool(to & threatenedByPawn) * 20000
                    : pt == ROOK
                      ? bool(to & threatenedByMinor) * 25000 + bool(to & threatenedByPawn) * 10000
                    : pt != PAWN ? bool(to & threatenedByPawn) * 15000
                                 : 0)
                 : 0),
            };

            for (int i = 0; i < N_HISTORY; ++i)
            {
                V += HISTORY_WEIGHT[i][0] * values[i] / HISTORY_WEIGHT[i][1];
                //                if (HISTORY_WEIGHT[i][0])
                //                   std::cerr << i << " => " << HISTORY_WEIGHT[i][0] * values[i] / HISTORY_WEIGHT[i][1]
                //                            << std::endl;
            }

            //std::cerr << values[HISTORY_EN_PRISE] << " " << V << std::endl;
            m.value += V;
            k++;
            if (MEASURE_BIAS)
            {
                int weight = getWeight(depth);
                dbg_mean_of(V, 0, weight);
                dbg_mean_of(V, depth, weight);
            }
        }
    }
}

// Returns the next move satisfying a predicate function.
// It never returns the TT move.
template<MovePicker::PickType T, typename Pred>
ExtMove MovePicker::select(Pred filter) {

    while (cur < endMoves)
    {
        if constexpr (T == Best)
            std::swap(*cur, *std::max_element(cur, endMoves));

        if (*cur != ttMove && filter())
            return *cur++;

        cur++;
    }
    return MOVE_NONE;
}

bool MovePicker::isQuiet() const { return stage == QUIET; }
bool MovePicker::isEvasion() const { return stage == EVASION; }
bool MovePicker::isRefutation() const { return stage == REFUTATION; }
bool MovePicker::isCapture() const { return stage == GOOD_CAPTURE || stage == BAD_CAPTURE; }

// Most important method of the MovePicker class. It
// returns a new pseudo-legal move every time it is called until there are no more
// moves left, picking the move with the highest score from a list of generated moves.
ExtMove MovePicker::next_move(bool skipQuiets) {

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
    case QCAPTURE_INIT :
        cur = endBadCaptures = moves;
        endMoves             = generate<CAPTURES>(pos, cur);

        score<SCORE_CAPTURES>();
        partial_insertion_sort(cur, endMoves, std::numeric_limits<int>::min());
        ++stage;
        goto top;

    case GOOD_CAPTURE :
        if (select<Next>([&]() {
                // Move losing capture to endBadCaptures to be tried later
                return pos.see_ge(*cur, Value(-cur->value)) ? true
                                                            : (*endBadCaptures++ = *cur, false);
            }))
            return *(cur - 1);

        // Prepare the pointers to loop over the refutations array
        cur      = std::begin(refutations);
        endMoves = std::end(refutations);

        // If the countermove is the same as a killer, skip it
        if (refutations[0].move == refutations[2].move
            || refutations[1].move == refutations[2].move)
            --endMoves;

        score<SCORE_REFUTATIONS>();

        ++stage;
        [[fallthrough]];

    case REFUTATION :
        if (select<Best>([&]() {
                return *cur != MOVE_NONE && !pos.capture_stage(*cur) && pos.pseudo_legal(*cur);
            }))
            return *(cur - 1);
        ++stage;
        [[fallthrough]];

    case QUIET_INIT :
        if (!skipQuiets)
        {
            cur      = endBadCaptures;
            endMoves = generate<QUIETS>(pos, cur);

            score<SCORE_QUIETS>();
            if (USE_FULL_QUIET_SORT)
                partial_insertion_sort(cur, endMoves, std::numeric_limits<int>::min());
            else
                partial_insertion_sort(cur, endMoves, -1960 - 3130 * depth /*+ 1114 + 24 * depth*/);
        }

        ++stage;
        [[fallthrough]];

    case QUIET :
        if (!skipQuiets && select<Next>([&]() {
                return *cur != refutations[0].move && *cur != refutations[1].move
                    && *cur != refutations[2].move;
            }))
            return *(cur - 1);

        // Prepare the pointers to loop over the bad captures
        cur      = moves;
        endMoves = endBadCaptures;

        ++stage;
        [[fallthrough]];

    case BAD_CAPTURE :
        return select<Next>([]() { return true; });

    case EVASION_INIT :
        cur      = moves;
        endMoves = generate<EVASIONS>(pos, cur);

        score<SCORE_EVASIONS>();
        ++stage;
        [[fallthrough]];

    case EVASION :
        return select<Best>([]() { return true; });

    case PROBCUT :
        return select<Next>([&]() { return pos.see_ge(*cur, threshold); });

    case QCAPTURE :
        if (select<Next>([]() { return true; }))
            return *(cur - 1);

        // If we did not find any move and we do not try checks, we have finished
        if (depth != DEPTH_QS_CHECKS)
            return MOVE_NONE;

        ++stage;
        [[fallthrough]];

    case QCHECK_INIT :
        cur      = moves;
        endMoves = generate<QUIET_CHECKS>(pos, cur);

        ++stage;
        [[fallthrough]];

    case QCHECK :
        return select<Next>([]() { return true; });
    }

    assert(false);
    return MOVE_NONE;  // Silence warning
}

}  // namespace Stockfish
