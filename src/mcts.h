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

#ifndef MCTS_H_INCLUDED
#define MCTS_H_INCLUDED

#include <map>

#include "types.h"

namespace Stockfish {

namespace MCTS {

class Tree;

class Node {

    static constexpr double C = 545;  // ~ 385 * sqrt(2)

    struct Stats {
        int64_t sum   = 0;
        int64_t count = 0;
        Node*   next  = nullptr;

        Value value() const { return count ? sum / count : VALUE_ZERO; }
    };

    Stats nodeStats;

    std::map<Move, Stats> moveStats;

   public:
    Node();
    void  clear();
    void  update(Move m, Value v);
    bool  hasValue(Move m) const;
    Value getValue(Move m) const;
    Value getUCB(Move m) const;
    Value getAverage() const;
    Node* next(Move m) const;

    friend Tree;
};

class Tree {
    static constexpr int SIZE = 16384;

    int                    nNodes = 1;
    std::array<Node, SIZE> nodes;

   public:
    Node*  root();
    void   clear();
    Node*  extendMove(Node* node, Move m);
    double usage() const;
};

}  // namespace MCTS

}  // namespace Stockfish

#endif  // #ifndef MCTS_H_INCLUDED
