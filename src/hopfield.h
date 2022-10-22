#ifndef HOPFIELD_H
#define HOPFIELD_H

#include <vector>
#include "types.h"

namespace Stockfish
{
    constexpr int BITS_PER_MOVE = 12;
    constexpr int MOVES_PER_PATTERN = 3;

    typedef std::vector<int> Pattern;

    class Hopfield 
    {
        static constexpr int RESOLUTION = 4096;

        int learnedPatterns;
        std::vector<std::vector<int>> weight;

        public:
        const int N;
        const int N_FIXED;

        Hopfield(int n, int n_fixed = 0);
        void clear();
        void addPattern(const Pattern& pattern);
        void addPattern(const std::vector<Pattern>& patterns);
        void retrievePattern(Pattern& input) const;
    };

    Hopfield::Hopfield(int n, int n_fixed) : learnedPatterns(0), weight(n, std::vector<int>(n, 0)), N(n), N_FIXED(n_fixed)
    {
    }

    void Hopfield::clear()
    {
        learnedPatterns = 0;
        for (auto w : weight)
            std::fill(w.begin(), w.end(), 0);
    }

    void Hopfield::addPattern(const Pattern& pattern)
    {
        ++learnedPatterns;

        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                if (i != j)
                   weight[i][j] = ((N - 1) * weight[i][j] + pattern[i] * pattern[j] * RESOLUTION) / N;

    }

    void Hopfield::addPattern(const std::vector<Pattern>& patterns)
    {
        for (auto p : patterns)
            addPattern(p);
    }

    void Hopfield::retrievePattern(Pattern& input) const
    {
        std::vector<int> index(N - N_FIXED);
        for (int i = 0; i < int(index.size()); ++i)
            index[i] = i;

        bool stable = false;
        while(!stable)
        {
            std::random_shuffle(index.begin(), index.end());

            stable = true;

            for(int i = 0; i < int(index.size()); ++i)
            {
                int sum = 0;
                for(int j = 0; j < N; ++j)
                    sum += input[j] * weight[index[i]][j];
                int out = (sum > 0 ? 1 : -1);
                if(out != input[index[i]])
                {
                    stable = false;
                    input[index[i]] = out;
                }
            }

        }
    }

    void buildPattern(Pattern& pattern, Move move, const std::vector<Move>& history)
    {
        for(int i = -1; i < int(history.size()); ++i)
        {
            int m = from_to((i < 0 ? move : history[i]));
            for(int j = 0; j < BITS_PER_MOVE; ++j)
                pattern[BITS_PER_MOVE * (i + 1) + j] = (m & (1 << j) ? 1 : -1);
        }
    }

    void setMove(Hopfield& hopfield, Move move, const std::vector<Move>& history)
    {
        Pattern pattern(hopfield.N);
        buildPattern(pattern, move, history);
        hopfield.addPattern(pattern);
    }

    Move extractMove(const Pattern& pattern, int moveNr = 0)
    {
        Move m = MOVE_NONE;
        for(int i = 0; i < BITS_PER_MOVE; ++i)
            m = Move(m * 2 + (pattern[BITS_PER_MOVE * moveNr + i] > 0 ? 1 : 0));
        return m;
    }

    Move getMove(const Hopfield& hopfield, const std::vector<Move>& history)
    {
        Pattern pattern(hopfield.N);
        buildPattern(pattern, MOVE_NONE, history);
        hopfield.retrievePattern(pattern);
        return extractMove(pattern);
    }

    Hopfield hopfield(BITS_PER_MOVE * MOVES_PER_PATTERN, BITS_PER_MOVE * (MOVES_PER_PATTERN - 1));
}

#endif
