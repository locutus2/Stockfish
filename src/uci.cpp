/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2022 The Stockfish developers (see AUTHORS file)

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

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <set>

#include "evaluate.h"
#include "movegen.h"
#include "position.h"
#include "search.h"
#include "thread.h"
#include "timeman.h"
#include "tt.h"
#include "uci.h"
#include "syzygy/tbprobe.h"

using namespace std;

namespace Stockfish {

extern vector<string> setup_bench(const Position&, istream&);
extern vector<string> setup_learn(const Position&, istream&, int& run);

namespace {

  // FEN string for the initial position in standard chess
  const char* StartFEN = "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1";


  // position() is called when the engine receives the "position" UCI command.
  // It sets up the position that is described in the given FEN string ("fen") or
  // the initial position ("startpos") and then makes the moves given in the following
  // move list ("moves").

  void position(Position& pos, istringstream& is, StateListPtr& states) {

    Move m;
    string token, fen;

    is >> token;

    if (token == "startpos")
    {
        fen = StartFEN;
        is >> token; // Consume the "moves" token, if any
    }
    else if (token == "fen")
        while (is >> token && token != "moves")
            fen += token + " ";
    else
        return;

    states = StateListPtr(new std::deque<StateInfo>(1)); // Drop the old state and create a new one
    pos.set(fen, Options["UCI_Chess960"], &states->back(), Threads.main());

    // Parse the move list, if any
    while (is >> token && (m = UCI::to_move(pos, token)) != MOVE_NONE)
    {
        states->emplace_back();
        pos.do_move(m, states->back());
    }
  }

  // trace_eval() prints the evaluation of the current position, consistent with
  // the UCI options set so far.

  void trace_eval(Position& pos) {

    StateListPtr states(new std::deque<StateInfo>(1));
    Position p;
    p.set(pos.fen(), Options["UCI_Chess960"], &states->back(), Threads.main());

    Eval::NNUE::verify();

    sync_cout << "\n" << Eval::trace(p) << sync_endl;
  }


  // setoption() is called when the engine receives the "setoption" UCI command.
  // The function updates the UCI option ("name") to the given value ("value").

  void setoption(istringstream& is) {

    string token, name, value;

    is >> token; // Consume the "name" token

    // Read the option name (can contain spaces)
    while (is >> token && token != "value")
        name += (name.empty() ? "" : " ") + token;

    // Read the option value (can contain spaces)
    while (is >> token)
        value += (value.empty() ? "" : " ") + token;

    if (Options.count(name))
        Options[name] = value;
    else
        sync_cout << "No such option: " << name << sync_endl;
  }


  // go() is called when the engine receives the "go" UCI command. The function
  // sets the thinking time and other parameters from the input string, then starts
  // with a search.

  void go(Position& pos, istringstream& is, StateListPtr& states) {

    Search::LimitsType limits;
    string token;
    bool ponderMode = false;

    limits.startTime = now(); // The search starts as early as possible

    while (is >> token)
        if (token == "searchmoves") // Needs to be the last command on the line
            while (is >> token)
                limits.searchmoves.push_back(UCI::to_move(pos, token));

        else if (token == "wtime")     is >> limits.time[WHITE];
        else if (token == "btime")     is >> limits.time[BLACK];
        else if (token == "winc")      is >> limits.inc[WHITE];
        else if (token == "binc")      is >> limits.inc[BLACK];
        else if (token == "movestogo") is >> limits.movestogo;
        else if (token == "depth")     is >> limits.depth;
        else if (token == "nodes")     is >> limits.nodes;
        else if (token == "movetime")  is >> limits.movetime;
        else if (token == "mate")      is >> limits.mate;
        else if (token == "perft")     is >> limits.perft;
        else if (token == "infinite")  limits.infinite = 1;
        else if (token == "ponder")    ponderMode = true;

    Threads.start_thinking(pos, states, limits, ponderMode);
  }


  // bench() is called when the engine receives the "bench" command.
  // Firstly, a list of UCI commands is set up according to the bench
  // parameters, then it is run one by one, printing a summary at the end.

  void bench(Position& pos, istream& args, StateListPtr& states) {

    string token;
    uint64_t num, nodes = 0, cnt = 1;

    vector<string> list = setup_bench(pos, args);
    num = count_if(list.begin(), list.end(), [](string s) { return s.find("go ") == 0 || s.find("eval") == 0; });

    TimePoint elapsed = now();

    for (const auto& cmd : list)
    {
        istringstream is(cmd);
        is >> skipws >> token;

        if (token == "go" || token == "eval")
        {
            cerr << "\nPosition: " << cnt++ << '/' << num << " (" << pos.fen() << ")" << endl;
            if (token == "go")
            {
               go(pos, is, states);
               Threads.main()->wait_for_search_finished();
               nodes += Threads.nodes_searched();
            }
            else
               trace_eval(pos);
        }
        else if (token == "setoption")  setoption(is);
        else if (token == "position")   position(pos, is, states);
        else if (token == "ucinewgame") { Search::clear(); elapsed = now(); } // Search::clear() may take a while
    }

    elapsed = now() - elapsed + 1; // Ensure positivity to avoid a 'divide by zero'

    dbg_print(); // Just before exiting
    dbg_printc(); // Just before exiting

    cerr << "\n==========================="
         << "\nTotal time (ms) : " << elapsed
         << "\nNodes searched  : " << nodes
         << "\nNodes/second    : " << 1000 * nodes / elapsed << endl;
  }

  void learn(Position& pos, istream& args, StateListPtr& states) {

    enum SCHEDULE {SCH_EXP, SCH_POLY, SCH_LIN};

    string token;
    uint64_t num, nodes = 0, cnt = 1;

    int run = 0;
    vector<string> list = setup_learn(pos, args, run);
    num = count_if(list.begin(), list.end(), [](string s) { return s.find("go ") == 0 || s.find("eval") == 0; });

    TimePoint elapsed = now();

    ACTIVE = false;
    for(int i = 0; i < N_PARAMS; ++i)
        PARAMS[i] = 0;

    std::vector<Move> bestMove;
    for (const auto& cmd : list)
    {
        istringstream is(cmd);
        is >> skipws >> token;

        if (token == "go" || token == "eval")
        {
            //cerr << "\nPosition: " << cnt++ << '/' << num << " (" << pos.fen() << ")" << endl;
            if (token == "go")
            {
               go(pos, is, states);
               Threads.main()->wait_for_search_finished();
               nodes += Threads.nodes_searched();
               bestMove.push_back(Threads.main()->rootMoves[0].pv[0]);
            }
            else
               trace_eval(pos);
        }
        else if (token == "setoption")  setoption(is);
        else if (token == "position")   position(pos, is, states);
        else if (token == "ucinewgame") { Search::clear(); elapsed = now(); } // Search::clear() may take a while
    }

    ACTIVE = true;
    //std::cerr << "nodes: " << nodes << std::endl;
    //std::cerr << "RUN: " << run << std::endl;
    std::srand(123456 + run);
    auto rngU = [](double a = 0, double b = 1)->double { return a + (b - a) * std::rand() / (double)RAND_MAX; };

    std::mt19937 gen{123456 + run};
    auto rngG = [&gen](double a = 0, double b = 1)->double { 
        std::normal_distribution<> d{(a+b)/2, (b-a)/4};
        return d(gen); 
    };

    auto energy = []()->double {
        return (  dbg_get_hit_on(0, 10, 10)
                - dbg_get_hit_on(1, 0, 10) + 1) / 2;
        //return dbg_get_hit_on();
    };

    constexpr double PCONT = 0.5;
    constexpr SCHEDULE schedule = SCH_EXP;
    constexpr bool DISCRETE = true;
    constexpr bool ONE_STEP = true;
    constexpr bool FULL_RANDOM = false;
    constexpr bool GAUSS = false;

    constexpr bool DYNAMIC_T0 = false;
    //constexpr double L = 10000;
    constexpr double POLY_ORDER = 10;
    constexpr double L = 0;
    constexpr double MIN_MOMENTUM = 0;
    constexpr double MAX_MOMENTUM = 0;//0.5;
    constexpr double DYNAMIC_MOM = 1;
    //constexpr double ALPHA = 0.001;
    //constexpr double ALPHA = 0.01; // base
    constexpr double ALPHA = 1;
    //constexpr double ALPHA = 1;
    //constexpr double T0 = 100000000;
    //constexpr double T0 = 2500000; // for ALPHA = 0.01
    //constexpr double T0 = 10000000; // for ALPHA = 0.1
    constexpr double ALPHA_BASE = 0.1;
    //constexpr double T_BASE = 10000000; // for ALPHA = 0.1
    double T_DIFF_MAX = 0;
    constexpr int KMAX = 1000;
    constexpr int LE = 1;//2 * N_PARAMS;//1;
    constexpr int RESTARTS = 1;
    double T_BASE = 1;//nodes; // for ALPHA = 0.1
    double T0 = T_BASE;// * std::pow(ALPHA / ALPHA_BASE, 0.6); // for ALPHA = 0.1
    //constexpr double BETA = POLY_TEMP ? POLY_ORDER : 0.98;
    //double BETA = schedule == SCH_POLY ? POLY_ORDER :
    //              schedule == SCH_LIN  ? (T0 - 1) / KMAX
    //          /* schedule == SCH_EXP */: std::pow(1 / T0, 1.0 / KMAX);
    constexpr double MIN_PARAM = -256;
    constexpr double MAX_PARAM = 256;
    constexpr int LOWER_PARAM = -1;
    constexpr int UPPER_PARAM = 1;
    //constexpr double MAX_PARAM = std::numeric_limits<double>::max();
    //constexpr int SHIFT = 128 * 4;
    //
    //double score0 = nodes;
    double score0 = energy();
    double score = score0;
    double POLD[N_PARAMS];
    double PBEST[N_PARAMS];
    double scorebest = score0;

    constexpr double P0 = 0.999;
    constexpr double DELTA0 = 1;
    T0 = -DELTA0 / std::log(P0);

    constexpr double P1 = 0.001;
    constexpr double T1 = -DELTA0 / std::log(P1);

    double BETA = schedule == SCH_POLY ? (1 - std::pow(T1/T0, 1.0/POLY_ORDER)) / KMAX:
                  schedule == SCH_LIN  ? (T0/T1 - 1) / KMAX
              /* schedule == SCH_EXP */: -std::log(T1 / T0) / KMAX;

    std::cerr << "BETA=" << BETA << std::endl;
    if (DYNAMIC_T0)
    {
        T0 = score0;
    }


    std::cerr << "=> BEST:";
    for(int i = 0; i < N_PARAMS; ++i)
    {
        std::cerr << " " << PARAMS[i];
    }
    std::cerr << std::endl;
    std::cerr << "it: " << 0 << " s: " << score << " T: " << T0 << " sr: " << score/score0 << std::endl;
    dbg_print();

    // init params

    for(int i = 0; i < N_PARAMS; ++i)
        PBEST[i] = PARAMS[i];

    for(int r =0; r <= RESTARTS; r++)
    {
        int it = 0;
        for(int k = 0; k < KMAX; k++)
        {
            //double T = schedule == SCH_POLY ? T0 * std::pow(1 - it / (double)KMAX, BETA) :
            double T = schedule == SCH_POLY ? T0 * std::pow(1 - BETA * k, POLY_ORDER) :
                       schedule == SCH_LIN  ? T0 / (1 + BETA * k)
                   /* schedule == SCH_EXP */: T0 * std::exp(-BETA * k); 
            for(int l = 0; l < LE; ++l)
            {
                dbg_clear();
                //std::cerr << "SET T=" << T << " T0=" << T0 << " T1=" << T1 << " it=" << it << std::endl;
                if (DISCRETE)
                {
                    for(int i = 0; i < N_PARAMS; ++i)
                        POLD[i] = PARAMS[i];

                    if (ONE_STEP)
                    {
                        int i = std::rand() % N_PARAMS;
                        PARAMS[i] = (PARAMS[i] - LOWER_PARAM + std::rand() % (UPPER_PARAM - LOWER_PARAM) + 1) 
                                   % (UPPER_PARAM - LOWER_PARAM + 1) + LOWER_PARAM;
                    }
                    else
                    {
                        std::set<int> used;
                        do
                        {
                            int i = std::rand() % N_PARAMS;
                            if(used.find(i) == used.end())
                            {
                                used.insert(i);
                                PARAMS[i] = (PARAMS[i] - LOWER_PARAM + std::rand() % (UPPER_PARAM - LOWER_PARAM) + 1)
                                           % (UPPER_PARAM - LOWER_PARAM + 1) + LOWER_PARAM;
                            }
                        }
                        while(rngU() <= PCONT);
                    }
                }
                else
                {
                    for(int i = 0; i < N_PARAMS; ++i)
                    {
                        POLD[i] = PARAMS[i];
                        double rnd;
                        if (FULL_RANDOM)
                        {
                            rnd = rngU(MIN_PARAM, MAX_PARAM);
                            PARAMS[i] = rnd;
                        }
                        else
                        {
                            if (GAUSS)
                                rnd = rngG(-1, 1);
                            else
                                rnd = rngU(-1, 1);

                            // random step to  a neighbour constraint by MIN_PARAM and MAX_PARAM
                            PARAMS[i] = std::min(MAX_PARAM, std::max(MIN_PARAM, PARAMS[i] + ALPHA * rnd));
                            // add also a momentum step in direction to the current best solution
                            PARAMS[i] += (MIN_MOMENTUM  + (MAX_MOMENTUM - MIN_MOMENTUM) * (1 - DYNAMIC_MOM * T / T0)) * (PBEST[i] - PARAMS[i]);
                            //PARAMS[i] += ALPHA * ((std::rand() % SHIFT) + (std::rand() % SHIFT ) - SHIFT + 1);
                        }
                    }
                }

                nodes = 0;
                std::vector<Move> bestMove2;
                for (const auto& cmd : list)
                {
                    istringstream is(cmd);
                    is >> skipws >> token;

                    if (token == "go" || token == "eval")
                    {
                        //cerr << "\nPosition: " << cnt++ << '/' << num << " (" << pos.fen() << ")" << endl;
                        if (token == "go")
                        {
                           go(pos, is, states);
                           Threads.main()->wait_for_search_finished();
                           nodes += Threads.nodes_searched();
                           bestMove2.push_back(Threads.main()->rootMoves[0].pv[0]);
                        }
                        else
                           trace_eval(pos);
                    }
                    else if (token == "setoption")  setoption(is);
                    else if (token == "position")   position(pos, is, states);
                    else if (token == "ucinewgame") { Search::clear(); elapsed = now(); } // Search::clear() may take a while
                }
                //std::cerr << "nodes: " << nodes << std::endl;

                //double new_score = nodes;
                double new_score = energy();
                for(int i = 0; i < (int)bestMove.size(); ++i)
                {
                    if (bestMove[i] != bestMove2[i])
                    {
                        new_score += L;
                    }
                }

                if(DYNAMIC_T0)
                {
                    //T_DIFF_MAX = std::max(T_DIFF_MAX, new_score - score);
                    T_DIFF_MAX = std::max(T_DIFF_MAX, std::abs(new_score - score));
                    T0 = 0.5 * T0 + 0.5 * T_DIFF_MAX;
                }

                //dbg_mean_of((new_score - score)/1024);
                //dbg_std_of((new_score - score)/1024);

                double delta = new_score - score;
                if (new_score < 1 && (delta < 0 || std::exp(-delta/T) >= rngU()))
                {
                    score = new_score;

                    if (score < scorebest)
                    {
                        scorebest = score;
                        for(int i = 0; i < N_PARAMS; ++i)
                            PBEST[i] = PARAMS[i];
                        std::cerr << "=> BEST";
                    }
                    else
                        std::cerr << "=> NEXT";
                    std::cerr << "[" << score/score0 << "]:";

                    for(int i = 0; i < N_PARAMS; ++i)
                    {
                            std::cerr << " " << PARAMS[i];
                    }
                    std::cerr << std::endl;
                }
                else
                {
                    for(int i = 0; i < N_PARAMS; ++i)
                    {
                        PARAMS[i] = POLD[i];
                    }
                }
                std::cerr << "it: " << it+1 << " s: " << score << " T: " << T << " sr: " << score/score0 << std::endl;
                dbg_print();
                it++;
            }
        }

        if (RESTARTS > 0)
        {
            std::cerr << "=> EPOCHE BEST[" << scorebest/score0 << "]:";
            for(int i = 0; i < N_PARAMS; ++i)
            {
                std::cerr << " " << PBEST[i];
            }
            std::cerr << std::endl;

            // set ne start point as best of last epoche
            for(int i = 0; i < N_PARAMS; ++i)
                PARAMS[i] = PBEST[i];
        }
    }

    std::cerr << "=> FINAL BEST[" << scorebest/score0 << "]:";
    for(int i = 0; i < N_PARAMS; ++i)
    {
        std::cerr << " " << PBEST[i];
    }
    std::cerr << std::endl;

    elapsed = now() - elapsed + 1; // Ensure positivity to avoid a 'divide by zero'

  }

  // The win rate model returns the probability of winning (in per mille units) given an
  // eval and a game ply. It fits the LTC fishtest statistics rather accurately.
  int win_rate_model(Value v, int ply) {

     // The model only captures up to 240 plies, so limit the input and then rescale
     double m = std::min(240, ply) / 64.0;

     // The coefficients of a third-order polynomial fit is based on the fishtest data
     // for two parameters that need to transform eval to the argument of a logistic
     // function.
     constexpr double as[] = {  -0.58270499,    2.68512549,   15.24638015,  344.49745382};
     constexpr double bs[] = {  -2.65734562,   15.96509799,  -20.69040836,   73.61029937 };

     // Enforce that NormalizeToPawnValue corresponds to a 50% win rate at ply 64
     static_assert(UCI::NormalizeToPawnValue == int(as[0] + as[1] + as[2] + as[3]));

     double a = (((as[0] * m + as[1]) * m + as[2]) * m) + as[3];
     double b = (((bs[0] * m + bs[1]) * m + bs[2]) * m) + bs[3];

     // Transform the eval to centipawns with limited range
     double x = std::clamp(double(v), -4000.0, 4000.0);

     // Return the win rate in per mille units rounded to the nearest value
     return int(0.5 + 1000 / (1 + std::exp((a - x) / b)));
  }

} // namespace


/// UCI::loop() waits for a command from the stdin, parses it and then calls the appropriate
/// function. It also intercepts an end-of-file (EOF) indication from the stdin to ensure a
/// graceful exit if the GUI dies unexpectedly. When called with some command-line arguments,
/// like running 'bench', the function returns immediately after the command is executed.
/// In addition to the UCI ones, some additional debug commands are also supported.

void UCI::loop(int argc, char* argv[]) {

  Position pos;
  string token, cmd;
  StateListPtr states(new std::deque<StateInfo>(1));

  pos.set(StartFEN, false, &states->back(), Threads.main());

  for (int i = 1; i < argc; ++i)
      cmd += std::string(argv[i]) + " ";

  do {
      if (argc == 1 && !getline(cin, cmd)) // Wait for an input or an end-of-file (EOF) indication
          cmd = "quit";

      istringstream is(cmd);

      token.clear(); // Avoid a stale if getline() returns nothing or a blank line
      is >> skipws >> token;

      if (    token == "quit"
          ||  token == "stop")
          Threads.stop = true;

      // The GUI sends 'ponderhit' to tell that the user has played the expected move.
      // So, 'ponderhit' is sent if pondering was done on the same move that the user
      // has played. The search should continue, but should also switch from pondering
      // to the normal search.
      else if (token == "ponderhit")
          Threads.main()->ponder = false; // Switch to the normal search

      else if (token == "uci")
          sync_cout << "id name " << engine_info(true)
                    << "\n"       << Options
                    << "\nuciok"  << sync_endl;

      else if (token == "setoption")  setoption(is);
      else if (token == "go")         go(pos, is, states);
      else if (token == "position")   position(pos, is, states);
      else if (token == "ucinewgame") Search::clear();
      else if (token == "isready")    sync_cout << "readyok" << sync_endl;

      // Add custom non-UCI commands, mainly for debugging purposes.
      // These commands must not be used during a search!
      else if (token == "flip")     pos.flip();
      else if (token == "bench")    bench(pos, is, states);
      else if (token == "learn")    learn(pos, is, states);
      else if (token == "d")        sync_cout << pos << sync_endl;
      else if (token == "eval")     trace_eval(pos);
      else if (token == "compiler") sync_cout << compiler_info() << sync_endl;
      else if (token == "export_net")
      {
          std::optional<std::string> filename;
          std::string f;
          if (is >> skipws >> f)
              filename = f;
          Eval::NNUE::save_eval(filename);
      }
      else if (token == "--help" || token == "help" || token == "--license" || token == "license")
          sync_cout << "\nStockfish is a powerful chess engine for playing and analyzing."
                       "\nIt is released as free software licensed under the GNU GPLv3 License."
                       "\nStockfish is normally used with a graphical user interface (GUI) and implements"
                       "\nthe Universal Chess Interface (UCI) protocol to communicate with a GUI, an API, etc."
                       "\nFor any further information, visit https://github.com/official-stockfish/Stockfish#readme"
                       "\nor read the corresponding README.md and Copying.txt files distributed along with this program.\n" << sync_endl;
      else if (!token.empty() && token[0] != '#')
          sync_cout << "Unknown command: '" << cmd << "'. Type help for more information." << sync_endl;

  } while (token != "quit" && argc == 1); // The command-line arguments are one-shot
}


/// UCI::value() converts a Value to a string by adhering to the UCI protocol specification:
///
/// cp <x>    The score from the engine's point of view in centipawns.
/// mate <y>  Mate in 'y' moves (not plies). If the engine is getting mated,
///           uses negative values for 'y'.

string UCI::value(Value v) {

  assert(-VALUE_INFINITE < v && v < VALUE_INFINITE);

  stringstream ss;

  if (abs(v) < VALUE_MATE_IN_MAX_PLY)
      ss << "cp " << v * 100 / NormalizeToPawnValue;
  else
      ss << "mate " << (v > 0 ? VALUE_MATE - v + 1 : -VALUE_MATE - v) / 2;

  return ss.str();
}


/// UCI::wdl() reports the win-draw-loss (WDL) statistics given an evaluation
/// and a game ply based on the data gathered for fishtest LTC games.

string UCI::wdl(Value v, int ply) {

  stringstream ss;

  int wdl_w = win_rate_model( v, ply);
  int wdl_l = win_rate_model(-v, ply);
  int wdl_d = 1000 - wdl_w - wdl_l;
  ss << " wdl " << wdl_w << " " << wdl_d << " " << wdl_l;

  return ss.str();
}


/// UCI::square() converts a Square to a string in algebraic notation (g1, a7, etc.)

std::string UCI::square(Square s) {
  return std::string{ char('a' + file_of(s)), char('1' + rank_of(s)) };
}


/// UCI::move() converts a Move to a string in coordinate notation (g1f3, a7a8q).
/// The only special case is castling where the e1g1 notation is printed in
/// standard chess mode and in e1h1 notation it is printed in Chess960 mode.
/// Internally, all castling moves are always encoded as 'king captures rook'.

string UCI::move(Move m, bool chess960) {

  Square from = from_sq(m);
  Square to = to_sq(m);

  if (m == MOVE_NONE)
      return "(none)";

  if (m == MOVE_NULL)
      return "0000";

  if (type_of(m) == CASTLING && !chess960)
      to = make_square(to > from ? FILE_G : FILE_C, rank_of(from));

  string move = UCI::square(from) + UCI::square(to);

  if (type_of(m) == PROMOTION)
      move += " pnbrqk"[promotion_type(m)];

  return move;
}


/// UCI::to_move() converts a string representing a move in coordinate notation
/// (g1f3, a7a8q) to the corresponding legal Move, if any.

Move UCI::to_move(const Position& pos, string& str) {

  if (str.length() == 5)
      str[4] = char(tolower(str[4])); // The promotion piece character must be lowercased

  for (const auto& m : MoveList<LEGAL>(pos))
      if (str == UCI::move(m, pos.is_chess960()))
          return m;

  return MOVE_NONE;
}

} // namespace Stockfish
