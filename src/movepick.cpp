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

std::vector<std::string> conditionName;

#define ADD_HISTORY_NEW(m,c)  ( \
	((m).history.size() >= conditionName.size() ? conditionName.push_back(#c) : void(0)), \
	(m).history.push_back((c)), \
	(c) \
)

#define ADD_HISTORY(m,c)  { \
	int hindex = (m).history.size(); \
	(m).history.push_back((c)); \
	if(hindex >= int(conditionName.size())) \
	    conditionName.push_back(#c); \
}

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

#ifdef USE_AVX512ICL
// Load the Move, and the ExtMove value, into all lanes of 512-bit registers
static void splat_extmove(const ExtMove& m, __m512i& move, __m512i& value) {
    move  = _mm512_set1_epi32(m.raw());
    value = _mm512_set1_epi32(m.value);
}

// Sorts up to 16 moves.
struct MoveSorter {
    static constexpr int MAX_ELEMENTS = 16;
    __m512i              sortedValues, sortedMoves;

    explicit MoveSorter(const ExtMove& first) {
        splat_extmove(first, sortedMoves, sortedValues);

        // Set the uninitialized move values to INT_MIN, so that they sort less than any other move
        sortedValues = _mm512_mask_set1_epi32(sortedValues, ~1, std::numeric_limits<int>::min());
    }

    void insert(const ExtMove& m) {
        __m512i move, value;
        splat_extmove(m, move, value);

        // Mask of all elements except the insertion point
        assert(m.value != std::numeric_limits<int>::min());
        const uint16_t expand = _kadd_mask16(_mm512_cmplt_epi32_mask(sortedValues, value), -1);

        sortedValues = _mm512_mask_expand_epi32(value, expand, sortedValues);
        sortedMoves  = _mm512_mask_expand_epi32(move, expand, sortedMoves);
    }

    void write_sorted(ExtMove* moves, std::ptrdiff_t count) const {
        static_assert(sizeof(ExtMove) == 8);
        assert(count <= MAX_ELEMENTS);

        // Because values and moves are stored separately, we need to reassemble the ExtMoves
        auto write = [&](int offset, const __m512i indices) {
            const __m512i extMoves = _mm512_permutex2var_epi32(sortedMoves, indices, sortedValues);
            const std::ptrdiff_t storeCount = count - offset;

            if (storeCount > 0)
                _mm512_mask_storeu_epi64(moves + offset, (1 << storeCount) - 1, extMoves);
        };

        write(0, _mm512_setr_epi32(0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23));
        write(8, _mm512_setr_epi32(8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15, 31));
    }
};
#endif

// Sort moves in descending order up to and including a given limit.
// The order of moves smaller than the limit is left unspecified.
void partial_insertion_sort(ExtMove* begin, ExtMove* end, int limit) {
    ExtMove *sortedEnd = begin, *p = begin + 1;

#ifdef USE_AVX512ICL
    if (begin == end)
        return;

    MoveSorter sorter(*begin);
    for (; p < end; ++p)
    {
        if (p->value >= limit)
        {
            if (sortedEnd - begin + 1 >= MoveSorter::MAX_ELEMENTS)  // sorter full
                break;

            sorter.insert(*p);
            *p = *++sortedEnd;
        }
    }
    sorter.write_sorted(begin, sortedEnd - begin + 1);
    // Use scalar implementation for any remaining elements
#endif

    for (; p < end; ++p)
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
                       const ButterflyHistory*      mh2,
                       const MaterialButterflyHistory*      mmh,
                       const LowPlyHistory*         lph,
                       const CapturePieceToHistory* cph,
                       const PieceToHistory**       ch,
                       const PieceToHistory**       seqh,
                       const SharedHistories*       sh,
                       int                          pl,
		       bool cc) :
    pos(p),
    mainHistory(mh),
    mainHistory2(mh2),
    materialMainHistory(mmh),
    lowPlyHistory(lph),
    captureHistory(cph),
    continuationHistory(ch),
    sequenceHistory(seqh),
    sharedHistory(sh),
    ttMove(ttm),
    depth(d),
    ply(pl),
    CC(cc){

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

constexpr bool MAX_THREATS = false;
constexpr bool MAX_CHECKS = false;
constexpr bool MAX_MAIN = false;

// Assigns a numerical value to each move in a list, used for sorting.
// Captures are ordered by Most Valuable Victim (MVV), preferring captures
// with a good history. Quiets moves are ordered using the history tables.
template<GenType Type>
ExtMove* MovePicker::score(const MoveList<Type>& ml) {

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
            //m.value = 2 * (*mainHistory)[us][m.raw()];
            m.value = 2 * (*mainHistory)[us][m.raw()] * (MAX_MAIN ? 59 : 32) / 32;
            m.value += 2 * sharedHistory->pawn_entry(pos)[pc][to];
            m.value += (*continuationHistory[0])[pc][to];
            m.value += (*continuationHistory[1])[pc][to];
            m.value += (*continuationHistory[2])[pc][to];
            m.value += (*continuationHistory[3])[pc][to];
            m.value += (*continuationHistory[5])[pc][to];

            // bonus for checks
            m.value += ((pos.check_squares(pt) & to) && pos.see_ge(m, -75)) * (MAX_CHECKS ? 40494 : 16384);
            //m.value += ((pos.check_squares(pt) & to) && pos.see_ge(m, -75)) * 16384;

            // penalty for moving to a square threatened by a lesser piece
            // or bonus for escaping an attack by a lesser piece.
            int v = (MAX_THREATS ? 4 : 20) * (bool(threatByLesser[pt] & from) - bool(threatByLesser[pt] & to));
            //int v = 20 * (bool(threatByLesser[pt] & from) - bool(threatByLesser[pt] & to));
            m.value += PieceValue[pt] * v;


            if (ply < LOW_PLY_HISTORY_SIZE)
                m.value += 8 * (*lowPlyHistory)[ply][m.raw()] / (1 + ply);

            m.history.clear();
	    if(CC)
	    {
		    ADD_HISTORY(m, 2*(*mainHistory)[us][m.raw()]);
		    ADD_HISTORY(m, 2*sharedHistory->pawn_entry(pos)[pc][to]);
		    ADD_HISTORY(m, (*continuationHistory[0])[pc][to]);
		    ADD_HISTORY(m, (*continuationHistory[1])[pc][to]);
		    ADD_HISTORY(m, (*continuationHistory[2])[pc][to]);
		    ADD_HISTORY(m, (*continuationHistory[3])[pc][to]);
		    ADD_HISTORY(m, (*continuationHistory[4])[pc][to]);
		    ADD_HISTORY(m, (*continuationHistory[5])[pc][to]);

		    ADD_HISTORY(m, ((pos.check_squares(pt) & to) && pos.see_ge(m, -75)) * 16384);
		    ADD_HISTORY(m, PieceValue[pt] * 20 * (bool(threatByLesser[pt] & from) - bool(threatByLesser[pt] & to)));
		    ADD_HISTORY(m, (ply < LOW_PLY_HISTORY_SIZE) * 8 * (*lowPlyHistory)[ply][m.raw()] / (1 + ply));

		    ADD_HISTORY(m, (*sequenceHistory[0])[pc][to]);
		    ADD_HISTORY(m, (*sequenceHistory[1])[pc][to]);
		    ADD_HISTORY(m, (*sequenceHistory[2])[pc][to]);
		    ADD_HISTORY(m, (*sequenceHistory[3])[pc][to]);
		    ADD_HISTORY(m, (*sequenceHistory[4])[pc][to]);
		    ADD_HISTORY(m, (*sequenceHistory[5])[pc][to]);

		    ADD_HISTORY(m, (pt == KNIGHT && more_than_one(Attacks::knight_attack(to) & pos.pieces(~us, KING, ROOK, QUEEN))) * 1024);
		    ADD_HISTORY(m, (pt == KNIGHT && bool(Attacks::knight_attack(to) & pos.pieces(~us, KING, ROOK, QUEEN))) * 1024);
		    ADD_HISTORY(m, (pt == KNIGHT) * (popcount(Attacks::knight_attack(to)) - popcount(Attacks::knight_attack(from))) * 1024);

		    ADD_HISTORY(m, (pt == PAWN && (to & (us == WHITE ? pawn_attacks_bb<BLACK>(pos.pieces(~us) ^ pos.pieces(~us, KING, PAWN))
					    : pawn_attacks_bb<WHITE>(pos.pieces(~us) ^ pos.pieces(~us, KING, PAWN))))) * 1024);

		    ADD_HISTORY(m, (*materialMainHistory)[us][m.raw()].get(material_index(pos)));
		    ADD_HISTORY(m, (*mainHistory2)[us][m.raw()]);
		    ADD_HISTORY(m, ((*mainHistory2)[us][m.raw()] + (*mainHistory)[us][m.raw()]) / 2);

		    int material = pos.count<ALL_PIECES>();
		    ADD_HISTORY(m, (*mainHistory)[us][m.raw()] * material / 32);
		    ADD_HISTORY(m, (*mainHistory)[us][m.raw()] * (32 - material) / 32);
		    ADD_HISTORY(m, sharedHistory->pawn_entry(pos)[pc][to] * material / 32);
		    ADD_HISTORY(m, sharedHistory->pawn_entry(pos)[pc][to] * (32 - material) / 32);

		    int piece_material = pos.count<ALL_PIECES>() - pos.count<PAWN>();
		    ADD_HISTORY(m, (*mainHistory)[us][m.raw()] * piece_material / 16);
		    ADD_HISTORY(m, (*mainHistory)[us][m.raw()] * (16 - piece_material) / 16);
		    ADD_HISTORY(m, sharedHistory->pawn_entry(pos)[pc][to] * piece_material / 16);
		    ADD_HISTORY(m, sharedHistory->pawn_entry(pos)[pc][to] * (16 - piece_material) / 16);

		    int pawn_material = pos.count<PAWN>();
		    ADD_HISTORY(m, (*mainHistory)[us][m.raw()] * pawn_material / 16);
		    ADD_HISTORY(m, (*mainHistory)[us][m.raw()] * (16 - pawn_material) / 16);
		    ADD_HISTORY(m, sharedHistory->pawn_entry(pos)[pc][to] * pawn_material / 16);
		    ADD_HISTORY(m, sharedHistory->pawn_entry(pos)[pc][to] * (16 - pawn_material) / 16);
	    }
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

    ExtMove e;
    e = Move::none();
    return e;
}

// This is the most important method of the MovePicker class. We emit one
// new pseudo-legal move on every call until there are no more moves left,
// picking the move with the highest score from a list of generated moves.
ExtMove MovePicker::next_move() {

    constexpr int compensation = - MAX_THREATS * 529 + MAX_CHECKS * 1610 - MAX_MAIN * 2535;
    constexpr int goodQuietThreshold = -14000 + compensation;
    ExtMove e;
top:
    switch (stage)
    {

    case MAIN_TT :
    case EVASION_TT :
    case QSEARCH_TT :
    case PROBCUT_TT :
        ++stage;
        e = ttMove;
        return e;

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

            partial_insertion_sort(cur, endCur, -3560 * depth + compensation);
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

        e = Move::none();
	return e;

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
        e = Move::none();
	return e;
}

void MovePicker::skip_quiet_moves() { skipQuiets = true; }

}  // namespace Stockfish
