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

    class Function;

    extern uint64_t START;

    class Term
    {
        static int nextId;

        protected:
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
        int generateId() const { return nextId++; }

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

        virtual Term* simplify() = 0;

        virtual ~Term()
        {
        }

        friend class Function;
    };

    template <int V>
    class Constant : public Term
    {
        static int id;

        public:
        Constant() : Term(id ? id : (id = generateId()), 0, std::to_string(V)) { }
        int operator()(const std::vector<int> &args) const { (void)args; return V; };
        Term* create() const { return new Constant<V>(); }
        Term* simplify() { return this; }
    };

    template <int I>
    class Variable : public Term
    {
        static int id;

        public:
        Variable() : Term(id ? id : (id = generateId()), 0, std::string("C") + std::to_string(I)) { }
        int operator()(const std::vector<int> &args) const { return args[I]; };
        Term* create() const { return new Variable<I>(); }
        Term* simplify() { return this; }
    };

    class UnaryFunction : public Term
    {
        protected:
        Term* operand;

        UnaryFunction(int id, const std::string& n, Term* op = nullptr) : Term(id, 1, n), operand(op) { }

        public:
        std::ostream& print(std::ostream& out) const
        {
            out << name;
            if (operand != nullptr)
            {
                out << (operand->getCount() > 1 ? "(" : "");
                operand->print(out);
                out << (operand->getCount() > 1 ? ")" : "");
            }
            return out;
        }

        void setOperand(Term* term, int i)
        {
            (void)i;
            operand = term;
        }

        Term* simplify()
        {
            operand = operand->simplify();
            UnaryFunction* next = dynamic_cast<UnaryFunction*>(operand);
            if (next != nullptr && code == next->getCode())
                return next->operand;
            return this;
        }
    };

    class Not : public UnaryFunction
    {
        static int id;

        public:
        Not(Term* op = nullptr) : UnaryFunction(id ? id : (id = generateId()), std::string("!"), op) { }
        int operator()(const std::vector<int> &args) const { return !(*operand)(args); };
        Term* create() const { return new Not(); }
    };

    class Negate : public UnaryFunction
    {
        static int id;

        public:
        Negate(Term* op = nullptr) : UnaryFunction(id ? id : (id = generateId()), std::string("-"), op) { }
        int operator()(const std::vector<int> &args) const { return -(*operand)(args); };
        Term* create() const { return new Negate(); }
    };

    class BinaryFunction : public Term
    {
        protected:
        Term* operand1;
        Term* operand2;

        BinaryFunction(int id, const std::string& n, Term* op1 = nullptr, Term* op2 = nullptr) : Term(id, 2, n), operand1(op1), operand2(op2) { }

        public:
        std::ostream& print(std::ostream& out) const
        {
            if (operand1 != nullptr && operand2 != nullptr)
            {
                bool parentOp1 = operand1->getCount() > 1 && code != operand1->getCode();
                bool parentOp2 = operand2->getCount() > 1 && code != operand2->getCode();
                out << (parentOp1 ? "(" : "");
                operand1->print(out);
                out << (parentOp1 ? ")" : "");
                out << name;
                out << (parentOp2 ? "(" : "");
                operand2->print(out);
                out << (parentOp2 ? ")" : "");
            }
            else
                out << name;

            return out;
        }

        void setOperand(Term* term, int i)
        {
            (i ? operand2 : operand1) = term;
        }

        Term* simplify()
        {
            operand1 = operand1->simplify();
            operand2 = operand2->simplify();
            return this;
        }
    };

    class And : public BinaryFunction
    {
        static int id;

        public:
        And(Term* op1 = nullptr, Term* op2 = nullptr) : BinaryFunction(id ? id : (id = generateId()), std::string("&&"), op1, op2) { }
        int operator()(const std::vector<int> &args) const { return (*operand1)(args) && (*operand2)(args); };
        Term* create() const { return new And(); }

        Term* simplify()
        {
            operand1 = operand1->simplify();
            operand2 = operand2->simplify();
            Constant<1>* next = dynamic_cast<Constant<1>*>(operand1);
            if (next != nullptr)
                return operand2;
            next = dynamic_cast<Constant<1>*>(operand2);
            if (next != nullptr)
                return operand1;
            return this;
        }
    };

    class Or : public BinaryFunction
    {
        static int id;

        public:
        Or(Term* op1 = nullptr, Term* op2 = nullptr) : BinaryFunction(id ? id : (id = generateId()), std::string("||"), op1, op2) { }
        int operator()(const std::vector<int> &args) const { return (*operand1)(args) || (*operand2)(args); };
        Term* create() const { return new Or(); }

        Term* simplify()
        {
            operand1 = operand1->simplify();
            operand2 = operand2->simplify();
            Constant<1>* next = dynamic_cast<Constant<1>*>(operand1);
            if (next != nullptr)
                return next;
            next = dynamic_cast<Constant<1>*>(operand2);
            if (next != nullptr)
                return next;
            return this;
        }
    };

    class Equal : public BinaryFunction
    {
        static int id;

        public:
        Equal(Term* op1 = nullptr, Term* op2 = nullptr) : BinaryFunction(id ? id : (id = generateId()), std::string("=="), op1, op2) { }
        int operator()(const std::vector<int> &args) const { return (*operand1)(args) == (*operand2)(args); };
        Term* create() const { return new Equal(); }
    };

    class Littler : public BinaryFunction
    {
        static int id;

        public:
        Littler(Term* op1 = nullptr, Term* op2 = nullptr) : BinaryFunction(id ? id : (id = generateId()), std::string("<"), op1, op2) { }
        int operator()(const std::vector<int> &args) const { return (*operand1)(args) < (*operand2)(args); };
        Term* create() const { return new Littler(); }
    };

    class Add : public BinaryFunction
    {
        static int id;

        public:
        Add(Term* op1 = nullptr, Term* op2 = nullptr) : BinaryFunction(id ? id : (id = generateId()), std::string("+"), op1, op2) { }
        int operator()(const std::vector<int> &args) const { return (*operand1)(args) + (*operand2)(args); };
        Term* create() const { return new Add(); }
    };

    class Subtract : public BinaryFunction
    {
        static int id;

        public:
        Subtract(Term* op1 = nullptr, Term* op2 = nullptr) : BinaryFunction(id ? id : (id = generateId()), std::string("-"), op1, op2) { }
        int operator()(const std::vector<int> &args) const { return (*operand1)(args) - (*operand2)(args); };
        Term* create() const { return new Subtract(); }
    };

    class Mult : public BinaryFunction
    {
        static int id;

        public:
        Mult(Term* op1 = nullptr, Term* op2 = nullptr) : BinaryFunction(id ? id : (id = generateId()), std::string("*"), op1, op2) { }
        int operator()(const std::vector<int> &args) const { return (*operand1)(args) * (*operand2)(args); };
        Term* create() const { return new Mult(); }

        Term* simplify()
        {
            operand1 = operand1->simplify();
            operand2 = operand2->simplify();
            Constant<1>* next = dynamic_cast<Constant<1>*>(operand1);
            if (next != nullptr)
                return operand2;
            next = dynamic_cast<Constant<1>*>(operand2);
            if (next != nullptr)
                return operand1;
            return this;
        }
    };

    class Div : public BinaryFunction
    {
        static int id;

        public:
        Div(Term* op1 = nullptr, Term* op2 = nullptr) : BinaryFunction(id ? id : (id = generateId()), std::string("/"), op1, op2) { }
        int operator()(const std::vector<int> &args) const { int denom = (*operand2)(args); return (denom ? (*operand1)(args) / denom : 0); };
        Term* create() const { return new Div(); }

        Term* simplify()
        {
            operand1 = operand1->simplify();
            operand2 = operand2->simplify();
            Constant<1>* next = dynamic_cast<Constant<1>*>(operand2);
            if (next != nullptr)
                return operand1;
            return this;
        }
    };

    class If : public Term
    {
        static int id;

        protected:
        Term* operand1;
        Term* operand2;
        Term* operand3;

        public:
        If(Term* op1 = nullptr, Term* op2 = nullptr, Term* op3 = nullptr) : Term(id ? id : (id = generateId()), 3, std::string("?")), operand1(op1), operand2(op2), operand3(op3) { }
        int operator()(const std::vector<int> &args) const { return (*operand1)(args) ? (*operand2)(args) : (*operand3)(args); };
        Term* create() const { return new If(); }

        std::ostream& print(std::ostream& out) const
        {
            if (operand1 != nullptr && operand2 != nullptr && operand3 != nullptr)
            {
                out << (operand1->getCount() > 1 ? "(" : "");
                operand1->print(out);
                out << (operand1->getCount() > 1 ? ")" : "");
                out << name;
                out << (operand2->getCount() > 1 ? "(" : "");
                operand2->print(out);
                out << (operand2->getCount() > 1 ? ")" : "");
                out << ':';
                out << (operand3->getCount() > 1 ? "(" : "");
                operand3->print(out);
                out << (operand3->getCount() > 1 ? ")" : "");
            }
            else
                out << name;

            return out;
        }


        void setOperand(Term* term, int i)
        {
            (i == 2 ? operand3 : i == 1 ? operand2 : operand1) = term;
        }

        Term* simplify()
        {
            operand1 = operand1->simplify();
            operand2 = operand2->simplify();
            operand3 = operand3->simplify();
            Constant<1>* next = dynamic_cast<Constant<1>*>(operand1);
            if (next != nullptr)
                return operand2;
            return this;
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

        void simplify()
        {
            root = root->simplify();
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
