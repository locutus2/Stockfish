#ifndef HOPFIELD_H
#define HOPFIELD_H

#include <vector>

#include "position.h"
#include "types.h"

namespace Stockfish
{
    class Hopfield
    {
        static constexpr int RESOLUTION = 4096;

        std::vector<std::vector<int>> weight;

        protected:
        typedef std::vector<int> Pattern;

        const int N;
        const int N_FIXED;

        public:

        Hopfield(int n, int n_fixed = 0) : weight(n, std::vector<int>(n, 0)), N(n), N_FIXED(n_fixed) {}
        void clear();
        void addPattern(const Pattern& pattern);
        void retrievePattern(Pattern& input) const;
    };

    class PositionHopfield : public Hopfield
    {
        static constexpr bool USE_POSITION = true;
        static constexpr bool ONE_HOT_ENCODING = true;
        static constexpr int BITS_PER_SQUARE = ONE_HOT_ENCODING ? 13 : 4;
        static constexpr int BITS_PER_POSITION = USE_POSITION ? 64 * BITS_PER_SQUARE : 0;
        static constexpr bool USE_INVERTED = false;
        static constexpr bool FIXED_HISTORY = true;
        static constexpr bool FULL_MOVE = false;
        static constexpr int BITS_PER_MOVE = 12 + 4 * FULL_MOVE;
        static constexpr int MOVES_PER_PATTERN = 1;

        public:
        PositionHopfield(int n = BITS_PER_MOVE * MOVES_PER_PATTERN + BITS_PER_POSITION, int n_fixed = (FIXED_HISTORY ? BITS_PER_MOVE * (MOVES_PER_PATTERN - 1) + BITS_PER_POSITION : 0)) : Hopfield(n, n_fixed) {}
        static void buildPattern(const Position& pos, Pattern& pattern, Move move, const std::vector<Move>& history);
        static Move extractMove(const Pattern& pattern, int moveNr = 0);
        void setMove(const Position& pos, Move move, const std::vector<Move>& history = std::vector<Move>());
        Move getMove(const Position& pos, const std::vector<Move>& history = std::vector<Move>()) const;
    };
}

#endif
