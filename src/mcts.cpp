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

#include "mcts.h"

namespace Stockfish {

namespace MCTS {

Node::Node() {}

void Node::update(Move m, Value v) {
    nodeStats.sum += v;
    nodeStats.count++;

    auto entry = moveStats.find(m);
    if (entry == moveStats.end())
        moveStats[m] = {v, 1};

    else
    {
        moveStats[m].sum += v;
        moveStats[m].count++;
    }
}

void Node::clear() {
    nodeStats.sum   = 0;
    nodeStats.count = 0;
    moveStats.clear();
}

Value Node::getValue(Move m) const {
    auto entry = moveStats.find(m);
    if (entry == moveStats.end())
        return VALUE_ZERO;
    else
        return entry->second.value();
}

Value Node::averageValue() const { return nodeStats.value(); }

Node* Node::next(Move m) const {
    auto entry = moveStats.find(m);
    if (entry == moveStats.end())
        return nullptr;
    else
        return entry->second.next;
}

Node* Tree::root() { return &nodes[0]; }

void Tree::clear() {
    nNodes = 1;
    nodes[0].clear();
};

Node* Tree::createNode() {
    if (nNodes < SIZE)
    {
        nodes[nNodes].clear();
        return &nodes[nNodes++];
    }
    else
        return nullptr;
}

double Tree::usage() const { return nNodes / double(SIZE); }

}  // namespace MCTS

}  // namespace Stockfish
