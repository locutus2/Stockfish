#include "hopfield.h"

namespace Stockfish
{
    void Hopfield::clear()
    {
        for (auto& w : weight)
            std::fill(w.begin(), w.end(), 0);
    }

    void Hopfield::addPattern(const Pattern& pattern)
    {
        const int WEIGHTS = N;
        for (int i = 0; i < N - N_FIXED; ++i)
            for (int j = i + 1; j < N; ++j)
                   weight[j][i] = weight[i][j] = ((WEIGHTS - 1) * weight[i][j] + pattern[i] * pattern[j] * RESOLUTION) / WEIGHTS;
    }

    void Hopfield::retrievePattern(Pattern& input) const
    {
        bool stable;
        std::vector<int> fixedSum(N - N_FIXED, 0);
        for(int i = 0; i < N - N_FIXED; ++i)
            for(int j = N - N_FIXED; j < N; ++j)
                fixedSum[i] += input[j] * weight[i][j];

        do
        {
            stable = true;

            for(int i = 0; i < N - N_FIXED; ++i)
            {
                int sum = fixedSum[i];
                for(int j = 0; j < N - N_FIXED; ++j)
                    sum += input[j] * weight[i][j];
                int out = (sum > 0 ? 1 : -1);
                if(out != input[i])
                {
                    stable = false;
                    input[i] = out;
                }
            }
        }
        while(!stable);
    }

    void PositionHopfield::buildPattern(const Position& pos, Pattern& pattern, Move move, const std::vector<Move>& history)
    {
        for(int i = -1; i < int(history.size()); ++i)
        {
            Move mm = (i < 0 ? move : history[i]);
            int m = FULL_MOVE ? mm : from_to(mm);
            for(int j = 0; j < BITS_PER_MOVE; ++j)
                pattern[BITS_PER_MOVE * (i + 1) + j] = (m & (1 << j) ? 1 : -1);
        }

        if (USE_POSITION)
        {
            constexpr int offset = BITS_PER_MOVE * MOVES_PER_PATTERN;

            if (ONE_HOT_ENCODING)
                std::fill(pattern.begin() + offset, pattern.end(),  -1);

            for (Square sq = SQ_A1; sq <= SQ_H8; ++sq)
            {
                int pieceIndex = pos.piece_on(sq);
                if (ONE_HOT_ENCODING)
                {
                    pieceIndex = pieceIndex - 2 * (pieceIndex > W_KING);
                    pattern[offset + BITS_PER_SQUARE * sq + pieceIndex] = 1;
                }
                else
                {
                    for (int i = 0; i < BITS_PER_SQUARE; ++i)
                        pattern[offset + BITS_PER_SQUARE * sq + i] = (pieceIndex & (1 << i) ? 1 : -1);

                }
            }
        }
    }

    Move PositionHopfield::extractMove(const Pattern& pattern, int moveNr)
    {
        Move m = MOVE_NONE;
        for(int i = 0; i < BITS_PER_MOVE; ++i)
            m = Move(m * 2 + (pattern[BITS_PER_MOVE * moveNr + i] > 0 ? 1 : 0));
        return m;
    }

    void PositionHopfield::setMove(const Position& pos, Move move, const std::vector<Move>& history)
    {
        assert(int(history.size()) == MOVES_PER_PATTERN - 1);

        Pattern pattern(N, 0);
        buildPattern(pos, pattern, move, history);
        addPattern(pattern);
    }

    Move PositionHopfield::getMove(const Position& pos, const std::vector<Move>& history) const
    {
        assert(int(history.size()) == MOVES_PER_PATTERN - 1);

        Pattern pattern(N, 0);
        buildPattern(pos, pattern, MOVE_NONE, history);
        retrievePattern(pattern);
        Move move = extractMove(pattern);

        if (USE_INVERTED)
        {
            if (is_ok(move) && pos.pseudo_legal(move))
                return move;

            else if (FULL_MOVE)
                move = Move(0xFFFF & ~move);

            else
                move = Move(from_to(Move(~move)));
        }

        return is_ok(move) ? move : MOVE_NONE;
    }
}
