/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2023 The Stockfish developers (see AUTHORS file)

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

#ifndef SEARCH_H_INCLUDED
#define SEARCH_H_INCLUDED

#include <string>
#include <vector>
#include <ostream>

#include "misc.h"
#include "movepick.h"
#include "types.h"

namespace Stockfish {

class Position;

namespace Learn {

    class Term
    {
        protected:
        static int id;

        const int code;
        const int count;
        const std::string name;

        public:
        Term(int c, int ac, const std::string& n) : code(c), count(ac), name(n)
        {
        }

        int getCode() { return code; }
        int getCount() { return count; }
        const std::string& getName() { return name; }
        int generateId() const { return id++; }

        virtual int operator()(const std::vector<int> &args) const = 0;
        virtual Term* create() const = 0;

        virtual std::ostream& print(std::ostream& out) const
        {
            return out << name;
        }

        virtual void setOperand(Term* term, int i)
        {
            (void)term;
            (void)i;
        }

        virtual ~Term()
        {
        }
    };

    template <int V>
    class Constant : public Term
    {
        static int id;

        public:
        Constant() : Term(id ? id : (id = generateId()), 0, std::string("v") + std::to_string(V))
        {
        }


        int operator()(const std::vector<int> &args) const { (void)args; return V; };

        Term* create() const
        {
            return new Constant<V>();
        }

        ~Constant()
        {
        }
    };

    template <int I>
    class Variable : public Term
    {
        static int id;

        public:
        Variable() : Term(id ? id : (id = generateId()), 0, std::string("x") + std::to_string(I))
        {
        }

        int operator()(const std::vector<int> &args) const { return args[I]; };

        Term* create() const
        {
            return new Variable<I>();
        }

        virtual ~Variable()
        {
        }
    };

    class Not : public Term
    {
        static int id;

        Term* operand;

        public:
        Not(Term* op = nullptr) : Term(id ? id : (id = generateId()), 1, std::string("!")), operand(op)
        {
        }

        std::ostream& print(std::ostream& out) const
        {
            out << name;
            if (operand != nullptr)
            {
                out << (operand->getCount() ? "(" : "");
                operand->print(out);
                out << (operand->getCount() ? ")" : "");
            }
            return out;
        }

        int operator()(const std::vector<int> &args) const { return !(*operand)(args); };

        void setOperand(Term* term, int i)
        {
            (void)i;
            operand = term;
        }

        Term* create() const
        {
            return new Not();
        }

        virtual ~Not()
        {
        }
    };

    class And : public Term
    {
        static int id;

        Term* operand1;
        Term* operand2;

        public:
        And(Term* op1 = nullptr, Term* op2 = nullptr) : Term(id ? id : (id = generateId()), 2, std::string("&&")), operand1(op1), operand2(op2)
        {
        }

        std::ostream& print(std::ostream& out) const
        {
            if (operand1 != nullptr && operand2 != nullptr)
            {
                out << (operand1->getCount() ? "(" : "");
                operand1->print(out);
                out << (operand1->getCount() ? ")" : "");
                out << name;
                out << (operand2->getCount() ? "(" : "");
                operand2->print(out);
                out << (operand2->getCount() ? ")" : "");
            }
            else
                out << name;

            return out;
        }

        int operator()(const std::vector<int> &args) const { return (*operand1)(args) && (*operand2)(args); };

        void setOperand(Term* term, int i)
        {
            (i ? operand2 : operand1) = term;
        }

        Term* create() const
        {
            return new And();
        }

        virtual ~And()
        {
        }
    };

    class Or : public Term
    {
        static int id;

        Term* operand1;
        Term* operand2;

        public:
        Or(Term* op1 = nullptr, Term* op2 = nullptr) : Term(id ? id : (id = generateId()), 2, std::string("||")), operand1(op1), operand2(op2)
        {
        }

        std::ostream& print(std::ostream& out) const
        {
            if (operand1 != nullptr && operand2 != nullptr)
            {
                out << (operand1->getCount() ? "(" : "");
                operand1->print(out);
                out << (operand1->getCount() ? ")" : "");
                out << name;
                out << (operand2->getCount() ? "(" : "");
                operand2->print(out);
                out << (operand2->getCount() ? ")" : "");
            }
            else
                out << name;

            return out;
        }

        int operator()(const std::vector<int> &args) const { return (*operand1)(args) || (*operand2)(args); };

        void setOperand(Term* term, int i)
        {
            (i ? operand2 : operand1) = term;
        }

        Term* create() const
        {
            return new Or();
        }

        virtual ~Or()
        {
        }
    };

    class Function
    {
        static bool initialized;
        static std::vector<Term*> functions;
        static int MaxF0;
        static int MaxF1;
        static int MaxF2;
        static int MaxF3;

        protected:
        Term *root;
         
        public:
        Function(uint64_t f)
        {
            if (!initialized)
            {
                initialized = true;
                initFunctions();
            }

            init(f);
        }

        bool operator()(const std::vector<int> &args) const { return (*root)(args); };

        std::ostream& print(std::ostream& out) const
        {
            return root->print(out);
        }

        static void initFunctions();
        void init(uint64_t f);

        
    };

    inline std::ostream& operator<<(std::ostream& out, const Function& f)
    {
        return f.print(out);
    }
}

namespace Search {


/// Stack struct keeps track of the information we need to remember from nodes
/// shallower and deeper in the tree during the search. Each search thread has
/// its own array of Stack objects, indexed by the current ply.

struct Stack {
  Move* pv;
  PieceToHistory* continuationHistory;
  int ply;
  Move currentMove;
  Move excludedMove;
  Move killers[2];
  Value staticEval;
  int statScore;
  int moveCount;
  bool inCheck;
  bool ttPv;
  bool ttHit;
  int doubleExtensions;
  int cutoffCnt;
};


/// RootMove struct is used for moves at the root of the tree. For each root move
/// we store a score and a PV (really a refutation in the case of moves which
/// fail low). Score is normally set at -VALUE_INFINITE for all non-pv moves.

struct RootMove {

  explicit RootMove(Move m) : pv(1, m) {}
  bool extract_ponder_from_tt(Position& pos);
  bool operator==(const Move& m) const { return pv[0] == m; }
  bool operator<(const RootMove& m) const { // Sort in descending order
    return m.score != score ? m.score < score
                            : m.previousScore < previousScore;
  }

  Value score = -VALUE_INFINITE;
  Value previousScore = -VALUE_INFINITE;
  Value averageScore = -VALUE_INFINITE;
  Value uciScore = -VALUE_INFINITE;
  bool scoreLowerbound = false;
  bool scoreUpperbound = false;
  int selDepth = 0;
  int tbRank = 0;
  Value tbScore;
  std::vector<Move> pv;
};

using RootMoves = std::vector<RootMove>;


/// LimitsType struct stores information sent by GUI about available time to
/// search the current move, maximum depth/time, or if we are in analysis mode.

struct LimitsType {

  LimitsType() { // Init explicitly due to broken value-initialization of non POD in MSVC
    time[WHITE] = time[BLACK] = inc[WHITE] = inc[BLACK] = npmsec = movetime = TimePoint(0);
    movestogo = depth = mate = perft = infinite = 0;
    nodes = 0;
  }

  bool use_time_management() const {
    return time[WHITE] || time[BLACK];
  }

  std::vector<Move> searchmoves;
  TimePoint time[COLOR_NB], inc[COLOR_NB], npmsec, movetime, startTime;
  int movestogo, depth, mate, perft, infinite;
  int64_t nodes;
};

extern LimitsType Limits;

void init();
void clear();

} // namespace Search

} // namespace Stockfish

#endif // #ifndef SEARCH_H_INCLUDED
