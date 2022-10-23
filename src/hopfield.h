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

        std::vector<std::vector<int>> weight;

        protected:
        const int N;
        const int N_FIXED;

        public:

        Hopfield(int n, int n_fixed = 0) : weight(n, std::vector<int>(n, 0)), N(n), N_FIXED(n_fixed) {}
        void clear();
        void addPattern(const Pattern& pattern);
        void retrievePattern(Pattern& input) const;
    };

    class MovesHopfield : public Hopfield
    {
        public:
        MovesHopfield(int n = BITS_PER_MOVE * MOVES_PER_PATTERN, int n_fixed = 0) : Hopfield(n, n_fixed) {}
        static void buildPattern(Pattern& pattern, Move move, const std::vector<Move>& history);
        static Move extractMove(const Pattern& pattern, int moveNr = 0);
        void setMove(Move move, const std::vector<Move>& history);
        Move getMove(const std::vector<Move>& history) const;
    };
}

#endif
