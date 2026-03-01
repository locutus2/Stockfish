/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2026 The Stockfish developers (see AUTHORS file)

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

bool MovePicker::isQuiet() const { return stage == GOOD_QUIET || stage == BAD_QUIET; }

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
                       const PieceToHistory*        ttmah,
                       const SharedHistories*       sh,
                       int                          pl,
                       Search::Stack*               s,
                       const std::vector<bool>&     cond,
                       NodeType                     nt) :
    pos(p),
    mainHistory(mh),
    lowPlyHistory(lph),
    captureHistory(cph),
    continuationHistory(ch),
    ttMoveAlternativeHistory(ttmah),
    sharedHistory(sh),
    ttMove(ttm),
    depth(d),
    ply(pl),
    ss(s),
    condition(cond),
    nodeType(nt) {

    if (pos.checkers())
        stage = EVASION_TT + !(ttm && pos.pseudo_legal(ttm));

    else
        stage = (depth > 0 ? MAIN_TT : QSEARCH_TT) + !(ttm && pos.pseudo_legal(ttm));
}

// MovePicker constructor for ProbCut: we generate captures with Static Exchange
// Evaluation (SEE) greater than or equal to the given threshold.
MovePicker::MovePicker(
  const Position& p, Move ttm, int th, const CapturePieceToHistory* cph, Search::Stack* s) :
    pos(p),
    captureHistory(cph),
    ttMove(ttm),
    threshold(th),
    ss(s) {
    assert(!pos.checkers());

    stage = PROBCUT_TT + !(ttm && pos.capture_stage(ttm) && pos.pseudo_legal(ttm));
}

constexpr int SCALE   = 128;
constexpr int W[2][7] = {
  //{-11,-3,17,-7,2,-10,8}, {11,-12,-7,-6,-2,-4,2}, // after 8 K
  //{-9,-8,9,-18,4,-13,5}, {11,-3,-7,-14,-4,-2,15}, // after 24 K
  //{-5, -2, 8, -15, -13, -15, 9}, {4, 5, -18, -19, -11, -15, -1},  // after 50 K
  //{-5,12,18,-29,-6,-22,13}, {2,8,-26,-30,-21,-11,-11}, // after 88 K
  {-8, 11, 14, -30, -5, -23, 5},
  {4, 4, -21, -33, -22, -6, -20},  // after 100 K
};

#define S(v) ((v) * 1)

constexpr int SCALE2    = 256;
constexpr int W2[10][3] = {
  //{ 4, 8, 7}, { -8, 8, -5}, { -2, 0, 1}, { -9, -5, 3}, { 6, 7, 1}, { 2, -1, 6}, { -1, 0, -5}, { 0, 6, 1}, { 1, 4, 6}, { -4, -5, 7} // after 35 K
  //{ 6, 5, 4}, { -7, 8, -7}, { -3, 0, 1}, { -7, -4, 0}, { 4, 2, -3}, { 3, -2, 7}, { -4, -1, -3}, { 2, 7, -1}, { 2, 3, 10}, { -4, -5, 6} // after 50 K
  //{ 9, 7, 0}, { -7, 7, -11}, { -6, -1, -1}, { -10, 1, 1}, { 4, 6, -2}, { 2, -4, 6}, { -7, -4, -3}, { 2, 10, -9}, { -2, 6, 9}, { 1, -8, 7} // after 87 K
  {9, 9, -2},   {-6, 11, -10}, {-7, 0, -1}, {-8, 1, 0}, {5, 6, -1}, {5, -4, 4},
  {-7, -4, -2}, {3, 11, -13},  {-2, 8, 10}, {1, -11, 8}  // after 100 K
  //{S(9), S(9), S(-2)},  {S(-6), S(11), S(-10)}, {S(-7), S(0), S(-1)},  {S(-8), S(1), S(0)}, {S(5), S(6), S(-1)},  {S(5), S(-4), S(4)},    {S(-7), S(-4), S(-2)}, {S(3), S(11), S(-13)}, {S(-2), S(8), S(10)}, {S(1), S(-11), S(8)}  // after 100 K double effect
};

// Assigns a numerical value to each move in a list, used for sorting.
// Captures are ordered by Most Valuable Victim (MVV), preferring captures
// with a good history. Quiets moves are ordered using the history tables.
template<GenType Type>
ExtMove* MovePicker::score(MoveList<Type>& ml) {

    static_assert(Type == CAPTURES || Type == QUIETS || Type == EVASIONS, "Wrong type");

    Color      us           = pos.side_to_move();
    const bool priorCapture = pos.captured_piece();
    const int  nt           = (nodeType == PV ? 1 : nodeType == CUT ? 2 : 0);

    [[maybe_unused]] Bitboard threatByLesser[KING + 1];
    if constexpr (Type == QUIETS)
    {
        threatByLesser[PAWN]   = 0;
        threatByLesser[KNIGHT] = threatByLesser[BISHOP] = pos.attacks_by<PAWN>(~us);
        threatByLesser[ROOK] =
          pos.attacks_by<KNIGHT>(~us) | pos.attacks_by<BISHOP>(~us) | threatByLesser[KNIGHT];
        threatByLesser[QUEEN] = pos.attacks_by<ROOK>(~us) | threatByLesser[ROOK];
        threatByLesser[KING]  = 0;
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
            m.value += 2 * sharedHistory->pawn_entry(pos)[pc][to];
            m.value += (*continuationHistory[0])[pc] [to];  // * (pos.captured_piece() && depth <= 1 ? 4 : 3) / 3;
            m.value += (*continuationHistory[1])[pc][to];
            m.value += (*continuationHistory[2])[pc][to];
            m.value += (*continuationHistory[3])[pc][to];
            m.value += (*continuationHistory[5])[pc][to];
            //m.value += (*continuationHistory[6])[pc][to];
            //m.value2 = (*ttMoveAlternativeHistory)[pc][to];
            //m.value2 = (*continuationHistory[0])[pc][to];
            //m.value2 = (*continuationHistory[1])[pc][to];
            //m.value2 = (*continuationHistory[2])[pc][to];
            //m.value2 = (*continuationHistory[3])[pc][to];
            //m.value2 = (*continuationHistory[5])[pc][to];
            //m.value2 = (*mainHistory)[us][m.raw()];
            //m.value2 = sharedHistory->pawn_entry(pos)[pc][to];

            /*
	    constexpr int S = 1;
	    constexpr int A = 1;
            m.value2 = 2 * (*mainHistory)[us][m.raw()] * (8 * S + A * std::max(7 - depth, 0)) / (8 * S);
            m.value2 += 2 * sharedHistory->pawn_entry(pos)[pc][to] * (8 * S + A * std::max(7 - depth, 0)) / (8 * S);
	   // * (8 * S + A * std::max(7 - depth, 0)) / (8 * S);
	   // * std::max(7 - depth, 1);
            m.value2 += (*continuationHistory[0])[pc][to];
            m.value2 += (*continuationHistory[1])[pc][to];// * (8 * S - A * std::max(7 - depth, 0)) / (8 * S);
            m.value2 += (*continuationHistory[2])[pc][to];
            m.value2 += (*continuationHistory[3])[pc][to];// * (8 * S - A * std::max(7 - depth, 0)) / (8 * S);
            m.value2 += (*continuationHistory[5])[pc][to];// * (8 * S - A * std::max(7 - depth, 0)) / (8 * S);
							  */
            /*
            m.value2 = 2 * (*mainHistory)[us][m.raw()] * (SCALE + W[priorCapture][0]);
            m.value2 += 2 * sharedHistory->pawn_entry(pos)[pc][to] * (SCALE + W[priorCapture][1]);
            m.value2 += (*continuationHistory[0])[pc][to] * (SCALE + W[priorCapture][2]);
            m.value2 += (*continuationHistory[1])[pc][to] * (SCALE + W[priorCapture][3]);
            m.value2 += (*continuationHistory[2])[pc][to] * (SCALE + W[priorCapture][4]);
            m.value2 += (*continuationHistory[3])[pc][to] * (SCALE + W[priorCapture][5]);
            m.value2 += (*continuationHistory[5])[pc][to] * (SCALE + W[priorCapture][6]);
            m.value2 /= SCALE;
	    */
            m.value2 = 2 * (*mainHistory)[us][m.raw()] * (SCALE2 + W2[0][nt]);
            m.value2 += 2 * sharedHistory->pawn_entry(pos)[pc][to] * (SCALE2 + W2[1][nt]);
            m.value2 += (*continuationHistory[0])[pc][to] * (SCALE2 + W2[2][nt]);
            m.value2 += (*continuationHistory[1])[pc][to] * (SCALE2 + W2[3][nt]);
            m.value2 += (*continuationHistory[2])[pc][to] * (SCALE2 + W2[4][nt]);
            m.value2 += (*continuationHistory[3])[pc][to] * (SCALE2 + W2[5][nt]);
            m.value2 += (*continuationHistory[5])[pc][to] * (SCALE2 + W2[6][nt]);

            // bonus for checks
            m.value += (bool(pos.check_squares(pt) & to) && pos.see_ge(m, -75)) * 16384;
            //m.value2 += (bool(pos.check_squares(pt) & to) && pos.see_ge(m, -75)) * 16384;
            m.value2 += (bool(pos.check_squares(pt) & to) && pos.see_ge(m, -75)) * 16384
                      * (SCALE2 + W[7][nt]);

            // penalty for moving to a square threatened by a lesser piece
            // or bonus for escaping an attack by a lesser piece.
            int v = 20 * (bool(threatByLesser[pt] & from) - bool(threatByLesser[pt] & to));
            m.value += PieceValue[pt] * v;
            //m.value2 += PieceValue[pt] * v;
            m.value2 += PieceValue[pt] * v * (SCALE2 + W2[8][nt]);


            if (ply < LOW_PLY_HISTORY_SIZE)
            {
                m.value += 8 * (*lowPlyHistory)[ply][m.raw()] / (1 + ply);
                //m.value2 += 8 * (*lowPlyHistory)[ply][m.raw()] / (1 + ply);
                m.value2 += 8 * (*lowPlyHistory)[ply][m.raw()] * (SCALE2 + W2[9][nt]) / (1 + ply);
            }
            m.value2 /= SCALE2;
            //m.value = m.value2;
            m.value2 = (*continuationHistory[6])[pc][to];
        }

        else  // Type == EVASIONS
        {
            if (pos.capture_stage(m))
                m.value = PieceValue[capturedPiece] + (1 << 28);
            else
                m.value = (*mainHistory)[us][m.raw()] + (*continuationHistory[0])[pc][to];
        }
    }
    return it;
}

// Returns the next move satisfying a predicate function.
// This never returns the TT move, as it was emitted before.
template<typename Pred>
ExtMove MovePicker::select(Pred filter) {

    for (; cur < endCur; ++cur)
        if (*cur != ttMove && filter())
            return *cur++;

    return Move::none();
}

// This is the most important method of the MovePicker class. We emit one
// new pseudo-legal move on every call until there are no more moves left,
// picking the move with the highest score from a list of generated moves.
ExtMove MovePicker::next_move() {

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
            //partial_insertion_sort(cur, endCur, std::numeric_limits<int>::min());
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
