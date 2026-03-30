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
#include <iostream>

#include "bitboard.h"
#include "misc.h"
#include "position.h"

namespace Stockfish {

constexpr double EV_LEN[5] = {2.1494756965236894717896241189081,
                           2.356734550000752630371545118127,
                           4.8427101039463429768825494126532,
                           1.8477751479486999529656430645896,
                           1.9317128633295891577058290228724};
constexpr double EV[5][6] = {
	{0.842415, 1.01132, 0.876589, 1.05802, 0, 1},
	{-1.94406, 0.0850495, -0.343454, 0.805999, 0, 1},
	{1.40054, 1.4072, -4.29948, 0.156808, 0, 1},
	{0.509304, -1.46188, -0.0759772, 0.109621, 0, 1},
	{-0.153756, 0.499881, 0.289916, -1.54076, 0, 1},
};

constexpr int SCALE_EV = 1024;
constexpr int SCALED_EV[5][6] = {
	{  401,  481,  417,  504, 0, 476 },
	{ -844,   36, -149,  350, 0, 434 },
	{  296,  297, -909,   33, 0, 211 },
	{  282, -810,  -42,   60, 0, 554 },
	{  -81,  264,  153, -816, 0, 530 }
};

#define S(i,j) (int(EV[(i)][(j)]/EV_LEN[(i)]*SCALE_EV))

/*
constexpr int SCALED_EV[5][6] = {
	{S(0,0),S(0,1),S(0,2),S(0,3),S(0,4),S(0,5)},
	{S(1,0),S(1,1),S(1,2),S(1,3),S(1,4),S(1,5)},
	{S(2,0),S(2,1),S(2,2),S(2,3),S(2,4),S(2,5)},
	{S(3,0),S(3,1),S(3,2),S(3,3),S(3,4),S(3,5)},
	{S(4,0),S(4,1),S(4,2),S(4,3),S(4,4),S(4,5)}
};
*/

constexpr int SCALE_W = 1024;
constexpr int CMH[6] = {SCALE_W,SCALE_W,SCALE_W,SCALE_W,0,SCALE_W}; // current master cmh weights

/*
 * CMH weights
3K 1012 1006 1029 1030 1037 
7K 1004 1006 1012 1038 1031
11K 1033 1003 979 1058 1023
17K 1023 975 996 1074 1037
20K 1036 957 1008 1063 1037
25K 1027 956 987 1076 1038
27K 1025 964 995 1079 1031
29K 1020 965 975 1080 1024
45K 1004 962 993 1046 1021
55K 995 989 1006 1025 991
77K 1004 991 1020 1005 989
83K 994 1006 1026 1014 992
100K 986 1001 1014 1030 1001
*/

/*
 * CMH weights C=25
 * 6K 1026 1030 1027 1020 1018
 * 11K 1016 1033 1027 1023 1019
 * 17K 1011 1034 1024 1026 1017
 * 30K 1016 1037 1024 1018 1019
 * 36K 1018 1033 1027 1016 1016
 * 50K 1016 1034 1026 1019 1020
 * 92K 1015 1031 1018 1024 1024
 * 100K 1012 1034 1017 1027 1021
 * */

/* CMH weights C=300
 * 21K 1079 964 973 1057 990
 * 46K 1042 995 1015 1088 930
 * 58K 1106 1024 1016 1113 973
 * 74K 1120 1024 970 1126 995
 * 96K 1109 1004 907 1147 1028
 */

//int W[5]={39,-10,128,54,-126};
int W[5];
int Random[5];

int w[6];

void init()
{
/*
	std::cerr << "SCALED_EV:" << std::endl;
	for(int i = 0; i < 5; i++)
	{
	     for(int j : {0,1,2,3,5})
		     std::cerr << SCALED_EV[i][j] << "\t";
	     std::cerr << std::endl;
 	}
*/
	//std::cerr << "CMH_WEIGHTS:" << std::endl;
        for(int c : {0,1,2,3,5})
	{
		w[c] = 0;
		for(int i = 0; i < 5; i++)
		    w[c] += W[i] * SCALED_EV[i][c];
		w[c] /= SCALE_EV;
		w[c] += CMH[c];
		//std::cerr << w[c] << " ";
	}
	//std::cerr << std::endl;
}

TUNE(SetRange(-SCALE_W,SCALE_W), W, Random, init);
UPDATE_ON_LAST();

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
                       const SharedHistories*       sh,
                       int                          pl) :
    pos(p),
    mainHistory(mh),
    lowPlyHistory(lph),
    captureHistory(cph),
    continuationHistory(ch),
    sharedHistory(sh),
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
            m.value = w[0] * (*continuationHistory[0])[pc][to];
            m.value += w[1] * (*continuationHistory[1])[pc][to];
            m.value += w[2] * (*continuationHistory[2])[pc][to];
            m.value += w[3] * (*continuationHistory[3])[pc][to];
            m.value += w[5] * (*continuationHistory[5])[pc][to];
            m.value /= SCALE_W;

            m.value += 2 * (*mainHistory)[us][m.raw()];
            m.value += 2 * sharedHistory->pawn_entry(pos)[pc][to];

            // bonus for checks
            m.value += (bool(pos.check_squares(pt) & to) && pos.see_ge(m, -75)) * 16384;

            // penalty for moving to a square threatened by a lesser piece
            // or bonus for escaping an attack by a lesser piece.
            int v = 20 * (bool(threatByLesser[pt] & from) - bool(threatByLesser[pt] & to));
            m.value += PieceValue[pt] * v;


            if (ply < LOW_PLY_HISTORY_SIZE)
                m.value += 8 * (*lowPlyHistory)[ply][m.raw()] / (1 + ply);
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
