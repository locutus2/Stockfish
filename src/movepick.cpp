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

int N_LGS;

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

// https://www.emathhelp.net/calculators/linear-algebra/eigenvalue-and-eigenvector-calculator
std::vector<double> pca(const std::vector<double>& X, const std::vector<std::vector<double>>& eVec)
{
	const int n = int(X.size());
	std::vector<double> v;
	for(int i = 0; i < int(eVec.size()); i++)
	{
		double s = 0;
	    for(int j = 0; j < n; j++)
		    s += X[j] * eVec[i][j];
	    v.push_back(s);
	}
	return v;
}

void lgs(const std::vector<double>& X, double Y, int offset = 0)
{
	const int n = int(X.size());
	N_LGS = n;
	for(int i = 0; i < n; i++)
	   for(int j = 0; j < n; j++)
	   {
		   dbg_sum_of(X[i] * X[j], offset+10*i+j);
	   }

	for(int i = 0; i < n; i++)
	{
	    dbg_sum_of(X[i] * Y, offset+100+i);
	}
}

/*
Correl. #0: Total 13279839 Coefficient 1
Correl. #1: Total 13279839 Coefficient 0.351996
Correl. #2: Total 13279839 Coefficient 0.0264105
Correl. #3: Total 13279839 Coefficient 0.0621023
Correl. #4: Total 13279839 Coefficient 0.0147296
Correl. #5: Total 13279839 Coefficient 0.0762523
Correl. #6: Total 13279839 Coefficient 0.0735225
Correl. #10: Total 13279839 Coefficient 0.351996
Correl. #11: Total 13279839 Coefficient 1
Correl. #12: Total 13279839 Coefficient 0.0610171
Correl. #13: Total 13279839 Coefficient 0.0912389
Correl. #14: Total 13279839 Coefficient 0.0500264
Correl. #15: Total 13279839 Coefficient 0.0998193
Correl. #16: Total 13279839 Coefficient 0.102667
Correl. #20: Total 13279839 Coefficient 0.0264105
Correl. #21: Total 13279839 Coefficient 0.0610171
Correl. #22: Total 13279839 Coefficient 1
Correl. #23: Total 13279839 Coefficient 0.286815
Correl. #24: Total 13279839 Coefficient 0.21815
Correl. #25: Total 13279839 Coefficient 0.268853
Correl. #26: Total 13279839 Coefficient 0.262502
Correl. #30: Total 13279839 Coefficient 0.0621023
Correl. #31: Total 13279839 Coefficient 0.0912389
Correl. #32: Total 13279839 Coefficient 0.286815
Correl. #33: Total 13279839 Coefficient 1
Correl. #34: Total 13279839 Coefficient 0.261919
Correl. #35: Total 13279839 Coefficient 0.339889
Correl. #36: Total 13279839 Coefficient 0.30987
Correl. #40: Total 13279839 Coefficient 0.0147296
Correl. #41: Total 13279839 Coefficient 0.0500264
Correl. #42: Total 13279839 Coefficient 0.21815
Correl. #43: Total 13279839 Coefficient 0.261919
Correl. #44: Total 13279839 Coefficient 1
Correl. #45: Total 13279839 Coefficient 0.333123
Correl. #46: Total 13279839 Coefficient 0.296777
Correl. #50: Total 13279839 Coefficient 0.0762523
Correl. #51: Total 13279839 Coefficient 0.0998193
Correl. #52: Total 13279839 Coefficient 0.268853
Correl. #53: Total 13279839 Coefficient 0.339889
Correl. #54: Total 13279839 Coefficient 0.333123
Correl. #55: Total 13279839 Coefficient 1
Correl. #56: Total 13279839 Coefficient 0.380507
Correl. #60: Total 13279839 Coefficient 0.0735225
Correl. #61: Total 13279839 Coefficient 0.102667
Correl. #62: Total 13279839 Coefficient 0.262502
Correl. #63: Total 13279839 Coefficient 0.30987
Correl. #64: Total 13279839 Coefficient 0.296777
Correl. #65: Total 13279839 Coefficient 0.380507
Correl. #66: Total 13279839 Coefficient 1

1
0.351996
0.0264105
0.0621023
0.0147296
0.0762523
0.0735225

0.351996
1
0.0610171
0.0912389
0.0500264
0.0998193
0.102667

0.0264105
0.0610171
1
0.286815
0.21815
0.268853
0.262502

0.0621023
0.0912389
0.286815
1
0.261919
0.339889
0.30987

0.0147296
0.0500264
0.21815
0.261919
1
0.333123
0.296777

0.0762523
0.0998193
0.268853
0.339889
0.333123
1
0.380507

0.0735225
0.102667
0.262502
0.30987
0.296777
0.380507
1

Wolfram alpha
{
	{1,         0.351996,  0.0264105, 0.0621023, 0.0147296, 0.0762523, 0.0735225},
        {0.351996,  1,         0.0610171, 0.0912389, 0.0500264, 0.0998193, 0.102667},
        {0.0264105, 0.0610171, 1,         0.286815,  0.21815,   0.268853,  0.262502},
	{0.0621023, 0.0912389, 0.286815, 1, 0.261919, 0.339889, 0.30987},
	{0.0147296, 0.0500264, 0.21815, 0.261919, 1, 0.333123, 0.296777},
	{0.0762523, 0.0998193, 0.268853, 0.339889, 0.333123, 1, 0.380507},
	{0.0735225, 0.102667, 0.262502, 0.30987, 0.296777, 0.380507, 1}
}

{ {1,         0.351996,  0.0264105, 0.0621023, 0.0147296, 0.0762523, 0.0735225}, {0.351996,  1,         0.0610171, 0.0912389, 0.0500264, 0.0998193,     0.102667}, {0.0264105, 0.0610171, 1,         0.286815,  0.21815,   0.268853,  0.262502}, {0.0621023, 0.0912389, 0.286815, 1, 0.261919, 0.339889, 0.30987}, {0.0147296, 0.0500264, 0.21815, 0.261919, 1, 0.333123, 0.296777}, {0.0762523, 0.0998193, 0.268853, 0.339889, 0.333123, 1, 0.380507}, {0.0735225,     0.102667, 0.262502, 0.30987, 0.296777, 0.380507, 1}}

corr of only cmh
{{1,0.286815,0.21815,0.268853,0.262502},{0.286815,1,0.261919,0.339889,0.30987},{0.21815,0.261919,1,0.333123,0.296777},{0.268853,0.339889,0.333123,1,0.380507},{0.262502,0.30987,0.296777,0.380507,1}}
 
ONly cmh hist
Input:
eigenvectors | (1 | 0.2868 | 0.218 | 0.2689 | 0.2625
0.2868 | 1 | 0.2619 | 0.3399 | 0.31
0.218 | 0.2619 | 1 | 0.3331 | 0.2968
0.2689 | 0.3399 | 0.3331 | 1 | 0.3805
0.2625 | 0.31 | 0.2968 | 0.3805 | 1)
Result:
v_1≈(0.855136, 0.964462, 0.909978, 1.04287, 1)
v_2≈(-4.38987, -1.35432, 3.05071, 1.23126, 1)
v_3≈(-1.40177, 1.40841, -1.9617, 0.599733, 1)
v_4≈(0.321395, -1.13675, -0.464457, 0.234128, 1)
v_5≈(-0.0785987, 0.385916, 0.363473, -1.5685, 1)
Corresponding eigenvalues:
λ_1≈2.19035
λ_2≈0.801766
λ_3≈0.71461
λ_4≈0.683208
λ_5≈0.610067	

--------------------
bench 16 1 16 pos1000.fen

Sum #0: Total 914374875 Average 2.40237e+08
Sum #100: Total 914374875 Average 2.50837e+08
Correl. #0: Total 914374875 Coefficient 1
Correl. #1: Total 914374875 Coefficient 0.132906
Correl. #2: Total 914374875 Coefficient 0.108003
Correl. #3: Total 914374875 Coefficient 0.115345
Correl. #4: Total 914374875 Coefficient 0.109248
Correl. #10: Total 914374875 Coefficient 0.132906
Correl. #11: Total 914374875 Coefficient 1
Correl. #12: Total 914374875 Coefficient 0.120423
Correl. #13: Total 914374875 Coefficient 0.171158
Correl. #14: Total 914374875 Coefficient 0.154068
Correl. #20: Total 914374875 Coefficient 0.108003
Correl. #21: Total 914374875 Coefficient 0.120423
Correl. #22: Total 914374875 Coefficient 1
Correl. #23: Total 914374875 Coefficient 0.140482
Correl. #24: Total 914374875 Coefficient 0.117659
Correl. #30: Total 914374875 Coefficient 0.115345
Correl. #31: Total 914374875 Coefficient 0.171158
Correl. #32: Total 914374875 Coefficient 0.140482
Correl. #33: Total 914374875 Coefficient 1
Correl. #34: Total 914374875 Coefficient 0.184791
Correl. #40: Total 914374875 Coefficient 0.109248
Correl. #41: Total 914374875 Coefficient 0.154068
Correl. #42: Total 914374875 Coefficient 0.117659
Correl. #43: Total 914374875 Coefficient 0.184791
Correl. #44: Total 914374875 Coefficient 1
Correl. #100: Total 914374875 Coefficient 1
Correl. #101: Total 914374875 Coefficient -0.0297756
Correl. #102: Total 914374875 Coefficient -0.0956135
Correl. #103: Total 914374875 Coefficient -0.0547608
Correl. #104: Total 914374875 Coefficient -0.00188724
Correl. #110: Total 914374875 Coefficient -0.0297756
Correl. #111: Total 914374875 Coefficient 1
Correl. #112: Total 914374875 Coefficient 0.00101075
Correl. #113: Total 914374875 Coefficient -0.0693329
Correl. #114: Total 914374875 Coefficient 0.0208709
Correl. #120: Total 914374875 Coefficient -0.0956135
Correl. #121: Total 914374875 Coefficient 0.00101075
Correl. #122: Total 914374875 Coefficient 1
Correl. #123: Total 914374875 Coefficient 0.0139638
Correl. #124: Total 914374875 Coefficient -0.0404246
Correl. #130: Total 914374875 Coefficient -0.0547608
Correl. #131: Total 914374875 Coefficient -0.0693329
Correl. #132: Total 914374875 Coefficient 0.0139638
Correl. #133: Total 914374875 Coefficient 1
Correl. #134: Total 914374875 Coefficient -0.0457782
Correl. #140: Total 914374875 Coefficient -0.00188724
Correl. #141: Total 914374875 Coefficient 0.0208709
Correl. #142: Total 914374875 Coefficient -0.0404246
Correl. #143: Total 914374875 Coefficient -0.0457782
Correl. #144: Total 914374875 Coefficient 1
Correl. #200: Total 914374875 Coefficient 1
Correl. #201: Total 914374875 Coefficient 0.998547
Correl. #210: Total 914374875 Coefficient 0.998547
Correl. #211: Total 914374875 Coefficient 1

===========================
Total time (ms) : 3712975
Nodes searched  : 182896450
Nodes/second    : 49258
I:
1
Inverse XtX:
4.16255e-09
XtY:
2.50837e+08
beta:
1.04412

beta:
-9.68742e-11    1.04268 -0.0140372      -0.0388051      -0.017137       0.0272964




Wolfram alpha: corr
{{1,0.132906,0.108003,0.115345,0.109248},{0.132906,1,0.120423,0.171158,0.154068},{0.108003,0.120423,1,0.140482,0.117659},{0.115345,0.171158,0.140482,1,0.184791},{0.109248,0.154068,0.117659,0.184791,1}}
reduced corr
{{1,0.13291,0.108,0.11535,0.10925},{0.13291,1,0.12042,0.17116,0.15407},{0.108,0.12042,1,0.14048,0.11766},{0.11535,0.17116,0.14048,1,0.18479},{0.10925,0.15407,0.11766,0.18479,1}}

eigenvectors | (1 | 0.13291 | 0.108 | 0.11535 | 0.10925
0.13291 | 1 | 0.12042 | 0.17116 | 0.15407
0.108 | 0.12042 | 1 | 0.14048 | 0.11766
0.11535 | 0.17116 | 0.14048 | 1 | 0.18479
0.10925 | 0.15407 | 0.11766 | 0.18479 | 1)
Result:
v_1≈(0.842415, 1.01132, 0.876589, 1.05802, 1)
v_2≈(-1.94406, 0.0850495, -0.343454, 0.805999, 1)
v_3≈(1.40054, 1.4072, -4.29948, 0.156808, 1)
v_4≈(0.509304, -1.46188, -0.0759772, 0.109621, 1)
v_5≈(-0.153756, 0.499881, 0.289916, -1.54076, 1)
Corresponding eigenvalues:
λ_1≈1.5465
λ_2≈0.909245
λ_3≈0.892917
λ_4≈0.841727
λ_5≈0.809613

beta(all):
1.47172e-10     1.03638 -0.0713816      -0.0142814      0.0237446       0.0255354

beta(pc1):
189.396 1.04387

	* */
void corr(const std::vector<double>& X, int offset = 0)
{
	const int n = int(X.size());
	for(int i = 0; i < n; i++)
	   for(int j = 0; j < n; j++)
	   {
		   dbg_correl_of(X[i], X[j], offset+10*i+j);
	   }
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
            m.value = (*continuationHistory[0])[pc][to];
            m.value += (*continuationHistory[1])[pc][to];
            m.value += (*continuationHistory[2])[pc][to];
            m.value += (*continuationHistory[3])[pc][to];
            m.value += (*continuationHistory[5])[pc][to];

	    std::vector<double> F = { 
		    (double)(*continuationHistory[0])[pc][to],
		    (double)(*continuationHistory[1])[pc][to],
		    (double)(*continuationHistory[2])[pc][to],
		    (double)(*continuationHistory[3])[pc][to],
		    (double)(*continuationHistory[5])[pc][to],
	    };
	    //corr(F);

	    std::vector<std::vector<double>> eVec = {
		    // bench 16 1 16 pos1000.fen
                    {0.842415, 1.01132, 0.876589, 1.05802, 1},
                    {-1.94406, 0.0850495, -0.343454, 0.805999, 1},
                    {1.40054, 1.4072, -4.29948, 0.156808, 1},
                    {0.509304, -1.46188, -0.0759772, 0.109621, 1},
                    {-0.153756, 0.499881, 0.289916, -1.54076, 1},
		    /* bench
		    {0.855136, 0.964462, 0.909978, 1.04287, 1},
		    {-4.38987, -1.35432, 3.05071, 1.23126, 1},
		    {-1.40177, 1.40841, -1.9617, 0.599733, 1},
		    {0.321395, -1.13675, -0.464457, 0.234128, 1},
		    {-0.0785987, 0.385916, 0.363473, -1.5685, 1},
		    */
		    //{1,0,0,0,0},
		    //{0,1,0,0,0},
		    //{0,0,1,0,0},
		    //{0,0,0,1,0},
		    //{0,0,0,0,1},
	    };
	    std::vector<double> X = pca(F, eVec);
	    //corr(X, 100);
	    double y = m.value;
	    std::vector<double> Y = {y, X[0]};
	    corr(Y, 200);

	    //X.resize(1);
	    X.insert(X.begin(), 1); // bias
	    lgs(X, y);

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
