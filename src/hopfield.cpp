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

    void MovesHopfield::buildPattern(Pattern& pattern, Move move, const std::vector<Move>& history)
    {
        for(int i = -1; i < int(history.size()); ++i)
        {
            Move mm = (i < 0 ? move : history[i]);
            int m = FULL_MOVE ? mm : from_to(mm);
            for(int j = 0; j < BITS_PER_MOVE; ++j)
                pattern[BITS_PER_MOVE * (i + 1) + j] = (m & (1 << j) ? 1 : -1);
        }
    }

    Move MovesHopfield::extractMove(const Pattern& pattern, int moveNr)
    {
        Move m = MOVE_NONE;
        for(int i = 0; i < BITS_PER_MOVE; ++i)
            m = Move(m * 2 + (pattern[BITS_PER_MOVE * moveNr + i] > 0 ? 1 : 0));
        return m;
    }

    void MovesHopfield::setMove(Move move, const std::vector<Move>& history)
    {
        Pattern pattern(N, 0);
        buildPattern(pattern, move, history);
        addPattern(pattern);
    }

    Move MovesHopfield::getMove(const Position& pos, const std::vector<Move>& history) const
    {
        Pattern pattern(N, 0);
        buildPattern(pattern, MOVE_NONE, history);
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
