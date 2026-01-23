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

#include <algorithm>
#include <cmath>

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

bool Node::hasValue(Move m) const {
    auto entry = moveStats.find(m);
    return entry != moveStats.end() && entry->second.count > 0;
}

Value Node::getValue(Move m) const {
    auto entry = moveStats.find(m);
    if (entry == moveStats.end())
        return VALUE_ZERO;
    else
        return entry->second.value();
}

Value Node::getUCB(Move m) const {
    auto entry = moveStats.find(m);
    if (nodeStats.count > 0 && entry != moveStats.end() && entry->second.count > 0)
        return double(entry->second.sum) / entry->second.count
             + C * std::sqrt(std::log(nodeStats.count) / entry->second.count);
    else
        return VALUE_MATE_IN_MAX_PLY;
}

Value Node::getAverage() const { return nodeStats.value(); }

Move Node::getBestMove() const {
    return moveStats.empty() ? Move::none()
                             : std::max_element(moveStats.begin(), moveStats.end(),
                                                [](const auto& s1, const auto& s2) {
                                                    return s1.second.value() < s2.second.value()
                                                        || (s1.second.value() == s2.second.value()
                                                            && s1.first < s2.first);
                                                })
                                 ->first;
}

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

Node* Tree::extendMove(Node* node, Move m) {
    if (nNodes < SIZE)
    {
        nodes[nNodes].clear();
        return node->moveStats[m].next = &nodes[nNodes++];
    }
    else
        return nullptr;
}

double Tree::usage() const { return nNodes / double(SIZE); }

}  // namespace MCTS

}  // namespace Stockfish
