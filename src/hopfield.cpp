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
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                if (i != j)
                   weight[i][j] = ((N - 1) * weight[i][j] + pattern[i] * pattern[j] * RESOLUTION) / N;

    }

    void Hopfield::retrievePattern(Pattern& input) const
    {
        bool stable = false;
        while(!stable)
        {
            stable = true;

            for(int i = 0; i < N - N_FIXED; ++i)
            {
                int sum = 0;
                for(int j = 0; j < N; ++j)
                    sum += input[j] * weight[i][j];
                int out = (sum > 0 ? 1 : -1);
                if(out != input[i])
                {
                    stable = false;
                    input[i] = out;
                }
            }

        }
    }

    void MovesHopfield::buildPattern(Pattern& pattern, Move move, const std::vector<Move>& history)
    {
        for(int i = -1; i < int(history.size()); ++i)
        {
            int m = from_to((i < 0 ? move : history[i]));
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

    Move MovesHopfield::getMove(const std::vector<Move>& history) const
    {
        Pattern pattern(N, 0);
        buildPattern(pattern, MOVE_NONE, history);
        retrievePattern(pattern);
        Move move = extractMove(pattern);
        return is_ok(move) ? move : MOVE_NONE;
    }
}
