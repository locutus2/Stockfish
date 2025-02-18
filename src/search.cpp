/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2025 The Stockfish developers (see AUTHORS file)

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

#include "search.h"

#include <algorithm>
#include <bitset>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <list>
#include <ratio>
#include <string>
#include <utility>
#include <vector>

#include "evaluate.h"
#include "history.h"
#include "misc.h"
#include "movegen.h"
#include "movepick.h"
#include "nnue/network.h"
#include "nnue/nnue_accumulator.h"
#include "nnue/nnue_common.h"
#include "nnue/nnue_misc.h"
#include "position.h"
#include "syzygy/tbprobe.h"
#include "thread.h"
#include "timeman.h"
#include "tt.h"
#include "uci.h"
#include "ucioption.h"

namespace Stockfish {

struct Expression {
    int expr    = 0;
    int good    = 0;
    int support = 0;

    inline double score() const { return support > 0 ? good / double(support) : 0.0; }
};

bool LEARNING_ENABLED = false;

std::vector<std::string> conditionNames = {
          "!PvNode",      
          "!ss->ttPv",      
          "capture",      
          "!capture",     
          "givesCheck",    
          "!givesCheck", 
          "ss->inCheck",
          "!ss->inCheck", 
          "priorCapture", 
          "!priorCapture", 
          "improving",   
          "!improving",
          "cutNode",
          "allNode",
          "ttCapture",
          "!ttCapture",
          "ttData.value > alpha",
          "ttData.value <= alpha",
          "correctionValue > -2322860",
          "correctionValue <= -2322860",
          "ttData.depth >= depth",
          "ttData.depth < depth",
          "ss->statScore > -5902",
          "ss->statScore <= -5902",
          "depth > 5",
          "depth <= 5",
          "moveCount > 6",
          "moveCount <= 6",
          "moveCount==1",
          "moveCount!=1",
          "move==ttData.move",
          "move!=ttData.move",
          "ttData.move",
          "!ttData.move",
          "ttHit",
          "!ttHit",
          "excludedMove",
          "!excludedMove",
          "(ss-1)->currentMove==Move::null()",
          "(ss-1)->currentMove!=Move::null()",
          "ttData.bound == BOUND_LOWER",
          "ttData.bound == BOUND_EXACT",
          "ttData.bound == BOUND_UPPER",
          "ttData.bound != BOUND_LOWER",
          "ttData.bound != BOUND_EXACT",
          "ttData.bound != BOUND_UPPER",
          "(ss+1)->cutoffCnt>3",
          "(ss+1)->cutoffCnt<=3",
          "type_of(movedPiece)==PAWN",
          "type_of(movedPiece)==KNIGHT",
          "type_of(movedPiece)==BISHOP",
          "type_of(movedPiece)==ROOK",
          "type_of(movedPiece)==QUEEN",
          "type_of(movedPiece)==KING",
          "type_of(movedPiece)!=PAWN",
          "type_of(movedPiece)!=KNIGHT",
          "type_of(movedPiece)!=BISHOP",
          "type_of(movedPiece)!=ROOK",
          "type_of(movedPiece)!=QUEEN",
          "type_of(movedPiece)!=KING",
          "!(ss-1)->ttPv",
          "(ss-1)->inCheck",
          "!(ss-1)->inCheck",
          "(ss-1)->ttHit",
          "!(ss-1)->ttHit",
          "(ss-1)->statScore > -5902",
          "(ss-1)->statScore <= -5902",
};

constexpr int NC = 10;  //18;

std::vector<int> conditionIndex;

std::vector<Expression> E(1 << NC), Enew;

void clearData()
{
    E = std::vector<Expression>(1 << NC);
}

void initSearchExpression() { 
    LEARNING_ENABLED = true; 

    for(int i = 0; i < int(conditionNames.size()); i++)
        conditionIndex.push_back(i);

    int seed = int(std::time(nullptr));
    //seed = 1739865092;
    std::cerr << "SEED: " << seed << std::endl;

    std::srand(seed);

    for(int i = int(conditionIndex.size()) - 1; i >= 0; i--)
    {
        int j = std::rand() % (i+1);
        if(i != j)
            std::swap(conditionIndex[i], conditionIndex[j]);
    }
}

void endSearchExpression() { LEARNING_ENABLED = false; }

void startIterationSearchExpression() { 
    dbg_clear();
    clearData();
}

void endIterationSearchExpression() { }

void searchExpression(int it, std::ostream& out) {
    out << "------------ Iteration " << it << " ---------------" << std::endl;

    for (int i = 0; i < int(E.size()); ++i)
    {
        E[i].expr = i;
    }

    const auto Eorig = E;

    std::sort(E.begin(), E.end(), [&](const auto& a, const auto& b) {
        return a.score() > b.score() || (a.score() == b.score() && a.support > b.support)
        || (a.score() == b.score() && a.support == b.support && popcount(a.expr) < popcount(b.expr))
        || (a.score() == b.score() && a.support == b.support && popcount(a.expr) == popcount(b.expr) && a.expr < b.expr);
    });

    constexpr bool DEBUG = false;

    if (DEBUG)
    {
        for (int i = 0; i < int(E.size()); ++i)
        {
            if (E[i].support > 0)
                out << "orig " << i << ". expr=" << std::bitset<NC>(E[i].expr)
                          << " score=" << 100 * E[i].score() << "% support=" << E[i].support
                          << " good=" << E[i].good << std::endl;
        }

        out << std::endl;
    }

    out << "Condition:" << std::endl;
    int worstValue = -1;
    int worstConditionIndex = -1;

    for (int k = NC - 1; k >= 0; --k)
    {
        for (int i = 0; i < int(E.size()); i++)
        {
            if (E[i].expr & (1 << k))
            {
                out << "Cond " << NC - 1 - k << " [" << conditionNames[conditionIndex[NC-1-k]] << "] : best " << i << ":"
                          << ". expr=" << std::bitset<NC>(E[i].expr)
                          << " score=" << 100 * E[i].score() << "% support=" << E[i].support
                          << " good=" << E[i].good << std::endl;
                if (i > worstValue)
                {
                    worstValue = i;
                    worstConditionIndex = NC-1-k;
                }
                break;
            }
        }
    }

    out << std::endl;
    out << "Best formulas:" << std::endl;

    for (int i = 0, best = -1, lastBest = -1; best < 0 && i < int(E.size()); i++)
    {
        if (lastBest < 0 || E[i].support > E[lastBest].support)
        {
            best = i;
        }

        if (best >= 0)
        {
            out << "best " << best << ":" << ". expr=" << std::bitset<NC>(E[best].expr)
                      << " score=" << 100 * E[best].score() << "% support=" << E[best].support
                      << " good=" << E[best].good << std::endl;
            lastBest = best;
            best     = -1;

            bool first = true;
            out << "=>";
            for (int k = NC - 1; k >= 0; --k)
                if(E[i].expr & (1 << k))
                {
                    if(!first) out << " &&";
                    first = false;
                    out << " (" << conditionNames[conditionIndex[NC-1-k]] << ")";
                }
            out << std::endl;
        }
    }

    if(worstConditionIndex >= 0)
    {
        out << "=> Replace condition " << worstConditionIndex << std::endl;
        conditionIndex.push_back(conditionIndex[worstConditionIndex]);
        conditionIndex.erase(conditionIndex.begin() + worstConditionIndex);
    }
}

int getIndex(const std::vector<bool>& C) {
    int index = 0;
    for (int i = 0; i < int(C.size()); i++)
        index = index * 2 + C[i];
    return index;
}

std::vector<bool>  filterConditions(const std::vector<bool>& C0) {
    std::vector<bool> C;
    for(int i = 0; i < NC; i++)
        C.push_back(C0[conditionIndex[i]]);
    return C;
}

void learn(bool T, const std::vector<bool>& C0) {
    if (!LEARNING_ENABLED)
        return;

    assert(C0.size() == conditionNames.size());

    std::vector<bool> C = filterConditions(C0);

    assert(C.size() == NC);

    int index = getIndex(C);

    assert(E.size() == (1 << NC));
    for (int i = 0; i < int(E.size()); i++)
    {
        assert(i < int(E.size()));
        assert(i >= 0);
        if ((index & i) == i)
        {
            E[i].support++;
            if (T)
                E[i].good++;
        }
    }
}

namespace TB = Tablebases;

void syzygy_extend_pv(const OptionsMap&            options,
                      const Search::LimitsType&    limits,
                      Stockfish::Position&         pos,
                      Stockfish::Search::RootMove& rootMove,
                      Value&                       v);

using namespace Search;

namespace {


// (*Scalers):
// The values with Scaler asterisks have proven non-linear scaling.
// They are optimized to time controls of 180 + 1.8 and longer,
// so changing them or adding conditions that are similar requires
// tests at these types of time controls.

// Futility margin
Value futility_margin(Depth d, bool noTtCutNode, bool improving, bool oppWorsening) {
    Value futilityMult       = 112 - 26 * noTtCutNode;
    Value improvingDeduction = improving * futilityMult * 2;
    Value worseningDeduction = oppWorsening * futilityMult / 3;

    return futilityMult * d - improvingDeduction - worseningDeduction;
}

constexpr int futility_move_count(bool improving, Depth depth) {
    return (3 + depth * depth) / (2 - improving);
}

int correction_value(const Worker& w, const Position& pos, const Stack* const ss) {
    const Color us    = pos.side_to_move();
    const auto  m     = (ss - 1)->currentMove;
    const auto  pcv   = w.pawnCorrectionHistory[pawn_structure_index<Correction>(pos)][us];
    const auto  micv  = w.minorPieceCorrectionHistory[minor_piece_index(pos)][us];
    const auto  wnpcv = w.nonPawnCorrectionHistory[WHITE][non_pawn_index<WHITE>(pos)][us];
    const auto  bnpcv = w.nonPawnCorrectionHistory[BLACK][non_pawn_index<BLACK>(pos)][us];
    const auto  cntcv =
      m.is_ok() ? (*(ss - 2)->continuationCorrectionHistory)[pos.piece_on(m.to_sq())][m.to_sq()]
                 : 0;

    return 6995 * pcv + 6593 * micv + 7753 * (wnpcv + bnpcv) + 6049 * cntcv;
}

// Add correctionHistory value to raw staticEval and guarantee evaluation
// does not hit the tablebase range.
Value to_corrected_static_eval(const Value v, const int cv) {
    return std::clamp(v + cv / 131072, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);
}

void update_correction_history(const Position& pos,
                               Stack* const    ss,
                               Search::Worker& workerThread,
                               const int       bonus) {
    const Move  m  = (ss - 1)->currentMove;
    const Color us = pos.side_to_move();

    static constexpr int nonPawnWeight = 165;

    workerThread.pawnCorrectionHistory[pawn_structure_index<Correction>(pos)][us]
      << bonus * 109 / 128;
    workerThread.minorPieceCorrectionHistory[minor_piece_index(pos)][us] << bonus * 141 / 128;
    workerThread.nonPawnCorrectionHistory[WHITE][non_pawn_index<WHITE>(pos)][us]
      << bonus * nonPawnWeight / 128;
    workerThread.nonPawnCorrectionHistory[BLACK][non_pawn_index<BLACK>(pos)][us]
      << bonus * nonPawnWeight / 128;

    if (m.is_ok())
        (*(ss - 2)->continuationCorrectionHistory)[pos.piece_on(m.to_sq())][m.to_sq()]
          << bonus * 138 / 128;
}

// History and stats update bonus, based on depth
int stat_bonus(Depth d) { return std::min(158 * d - 98, 1622); }

// History and stats update malus, based on depth
int stat_malus(Depth d) { return std::min(802 * d - 243, 2850); }

// Add a small random component to draw evaluations to avoid 3-fold blindness
Value value_draw(size_t nodes) { return VALUE_DRAW - 1 + Value(nodes & 0x2); }
Value value_to_tt(Value v, int ply);
Value value_from_tt(Value v, int ply, int r50c);
void  update_pv(Move* pv, Move move, const Move* childPv);
void  update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus);
void  update_quiet_histories(
   const Position& pos, Stack* ss, Search::Worker& workerThread, Move move, int bonus);
void update_all_stats(const Position&      pos,
                      Stack*               ss,
                      Search::Worker&      workerThread,
                      Move                 bestMove,
                      Square               prevSq,
                      ValueList<Move, 32>& quietsSearched,
                      ValueList<Move, 32>& capturesSearched,
                      Depth                depth,
                      bool                 isTTMove,
                      int                  moveCount);

}  // namespace

Search::Worker::Worker(SharedState&                    sharedState,
                       std::unique_ptr<ISearchManager> sm,
                       size_t                          threadId,
                       NumaReplicatedAccessToken       token) :
    // Unpack the SharedState struct into member variables
    threadIdx(threadId),
    numaAccessToken(token),
    manager(std::move(sm)),
    options(sharedState.options),
    threads(sharedState.threads),
    tt(sharedState.tt),
    networks(sharedState.networks),
    refreshTable(networks[token]) {
    clear();
}

void Search::Worker::ensure_network_replicated() {
    // Access once to force lazy initialization.
    // We do this because we want to avoid initialization during search.
    (void) (networks[numaAccessToken]);
}

void Search::Worker::start_searching() {

    // Non-main threads go directly to iterative_deepening()
    if (!is_mainthread())
    {
        iterative_deepening();
        return;
    }

    main_manager()->tm.init(limits, rootPos.side_to_move(), rootPos.game_ply(), options,
                            main_manager()->originalTimeAdjust);
    tt.new_search();

    if (rootMoves.empty())
    {
        rootMoves.emplace_back(Move::none());
        main_manager()->updates.onUpdateNoMoves(
          {0, {rootPos.checkers() ? -VALUE_MATE : VALUE_DRAW, rootPos}});
    }
    else
    {
        threads.start_searching();  // start non-main threads
        iterative_deepening();      // main thread start searching
    }

    // When we reach the maximum depth, we can arrive here without a raise of
    // threads.stop. However, if we are pondering or in an infinite search,
    // the UCI protocol states that we shouldn't print the best move before the
    // GUI sends a "stop" or "ponderhit" command. We therefore simply wait here
    // until the GUI sends one of those commands.
    while (!threads.stop && (main_manager()->ponder || limits.infinite))
    {}  // Busy wait for a stop or a ponder reset

    // Stop the threads if not already stopped (also raise the stop if
    // "ponderhit" just reset threads.ponder)
    threads.stop = true;

    // Wait until all threads have finished
    threads.wait_for_search_finished();

    // When playing in 'nodes as time' mode, subtract the searched nodes from
    // the available ones before exiting.
    if (limits.npmsec)
        main_manager()->tm.advance_nodes_time(threads.nodes_searched()
                                              - limits.inc[rootPos.side_to_move()]);

    Worker* bestThread = this;
    Skill   skill =
      Skill(options["Skill Level"], options["UCI_LimitStrength"] ? int(options["UCI_Elo"]) : 0);

    if (int(options["MultiPV"]) == 1 && !limits.depth && !limits.mate && !skill.enabled()
        && rootMoves[0].pv[0] != Move::none())
        bestThread = threads.get_best_thread()->worker.get();

    main_manager()->bestPreviousScore        = bestThread->rootMoves[0].score;
    main_manager()->bestPreviousAverageScore = bestThread->rootMoves[0].averageScore;

    // Send again PV info if we have a new best thread
    if (bestThread != this)
        main_manager()->pv(*bestThread, threads, tt, bestThread->completedDepth);

    std::string ponder;

    if (bestThread->rootMoves[0].pv.size() > 1
        || bestThread->rootMoves[0].extract_ponder_from_tt(tt, rootPos))
        ponder = UCIEngine::move(bestThread->rootMoves[0].pv[1], rootPos.is_chess960());

    auto bestmove = UCIEngine::move(bestThread->rootMoves[0].pv[0], rootPos.is_chess960());
    main_manager()->updates.onBestmove(bestmove, ponder);
}

// Main iterative deepening loop. It calls search()
// repeatedly with increasing depth until the allocated thinking time has been
// consumed, the user stops the search, or the maximum search depth is reached.
void Search::Worker::iterative_deepening() {

    SearchManager* mainThread = (is_mainthread() ? main_manager() : nullptr);

    Move pv[MAX_PLY + 1];

    Depth lastBestMoveDepth = 0;
    Value lastBestScore     = -VALUE_INFINITE;
    auto  lastBestPV        = std::vector{Move::none()};

    Value  alpha, beta;
    Value  bestValue     = -VALUE_INFINITE;
    Color  us            = rootPos.side_to_move();
    double timeReduction = 1, totBestMoveChanges = 0;
    int    delta, iterIdx                        = 0;

    // Allocate stack with extra size to allow access from (ss - 7) to (ss + 2):
    // (ss - 7) is needed for update_continuation_histories(ss - 1) which accesses (ss - 6),
    // (ss + 2) is needed for initialization of cutOffCnt.
    Stack  stack[MAX_PLY + 10] = {};
    Stack* ss                  = stack + 7;

    for (int i = 7; i > 0; --i)
    {
        (ss - i)->continuationHistory =
          &this->continuationHistory[0][0][NO_PIECE][0];  // Use as a sentinel
        (ss - i)->continuationCorrectionHistory = &this->continuationCorrectionHistory[NO_PIECE][0];
        (ss - i)->staticEval                    = VALUE_NONE;
        (ss - i)->reduction                     = 0;
    }

    for (int i = 0; i <= MAX_PLY + 2; ++i)
    {
        (ss + i)->ply       = i;
        (ss + i)->reduction = 0;
    }

    ss->pv = pv;

    if (mainThread)
    {
        if (mainThread->bestPreviousScore == VALUE_INFINITE)
            mainThread->iterValue.fill(VALUE_ZERO);
        else
            mainThread->iterValue.fill(mainThread->bestPreviousScore);
    }

    size_t multiPV = size_t(options["MultiPV"]);
    Skill skill(options["Skill Level"], options["UCI_LimitStrength"] ? int(options["UCI_Elo"]) : 0);

    // When playing with strength handicap enable MultiPV search that we will
    // use behind-the-scenes to retrieve a set of possible moves.
    if (skill.enabled())
        multiPV = std::max(multiPV, size_t(4));

    multiPV = std::min(multiPV, rootMoves.size());

    int searchAgainCounter = 0;

    lowPlyHistory.fill(95);

    // Iterative deepening loop until requested to stop or the target depth is reached
    while (++rootDepth < MAX_PLY && !threads.stop
           && !(limits.depth && mainThread && rootDepth > limits.depth))
    {
        // Age out PV variability metric
        if (mainThread)
            totBestMoveChanges /= 2;

        // Save the last iteration's scores before the first PV line is searched and
        // all the move scores except the (new) PV are set to -VALUE_INFINITE.
        for (RootMove& rm : rootMoves)
            rm.previousScore = rm.score;

        size_t pvFirst = 0;
        pvLast         = 0;

        if (!threads.increaseDepth)
            searchAgainCounter++;

        // MultiPV loop. We perform a full root search for each PV line
        for (pvIdx = 0; pvIdx < multiPV; ++pvIdx)
        {
            if (pvIdx == pvLast)
            {
                pvFirst = pvLast;
                for (pvLast++; pvLast < rootMoves.size(); pvLast++)
                    if (rootMoves[pvLast].tbRank != rootMoves[pvFirst].tbRank)
                        break;
            }

            // Reset UCI info selDepth for each depth and each PV line
            selDepth = 0;

            // Reset aspiration window starting size
            delta     = 5 + std::abs(rootMoves[pvIdx].meanSquaredScore) / 13000;
            Value avg = rootMoves[pvIdx].averageScore;
            alpha     = std::max(avg - delta, -VALUE_INFINITE);
            beta      = std::min(avg + delta, VALUE_INFINITE);

            // Adjust optimism based on root move's averageScore
            optimism[us]  = 138 * avg / (std::abs(avg) + 81);
            optimism[~us] = -optimism[us];

            // Start with a small aspiration window and, in the case of a fail
            // high/low, re-search with a bigger window until we don't fail
            // high/low anymore.
            int failedHighCnt = 0;
            while (true)
            {
                // Adjust the effective depth searched, but ensure at least one
                // effective increment for every four searchAgain steps (see issue #2717).
                Depth adjustedDepth =
                  std::max(1, rootDepth - failedHighCnt - 3 * (searchAgainCounter + 1) / 4);
                rootDelta = beta - alpha;
                bestValue = search<Root>(rootPos, ss, alpha, beta, adjustedDepth, false);

                // Bring the best move to the front. It is critical that sorting
                // is done with a stable algorithm because all the values but the
                // first and eventually the new best one is set to -VALUE_INFINITE
                // and we want to keep the same order for all the moves except the
                // new PV that goes to the front. Note that in the case of MultiPV
                // search the already searched PV lines are preserved.
                std::stable_sort(rootMoves.begin() + pvIdx, rootMoves.begin() + pvLast);

                // If search has been stopped, we break immediately. Sorting is
                // safe because RootMoves is still valid, although it refers to
                // the previous iteration.
                if (threads.stop)
                    break;

                // When failing high/low give some update before a re-search. To avoid
                // excessive output that could hang GUIs like Fritz 19, only start
                // at nodes > 10M (rather than depth N, which can be reached quickly)
                if (mainThread && multiPV == 1 && (bestValue <= alpha || bestValue >= beta)
                    && nodes > 10000000)
                    main_manager()->pv(*this, threads, tt, rootDepth);

                // In case of failing low/high increase aspiration window and re-search,
                // otherwise exit the loop.
                if (bestValue <= alpha)
                {
                    beta  = (alpha + beta) / 2;
                    alpha = std::max(bestValue - delta, -VALUE_INFINITE);

                    failedHighCnt = 0;
                    if (mainThread)
                        mainThread->stopOnPonderhit = false;
                }
                else if (bestValue >= beta)
                {
                    beta = std::min(bestValue + delta, VALUE_INFINITE);
                    ++failedHighCnt;
                }
                else
                    break;

                delta += delta / 3;

                assert(alpha >= -VALUE_INFINITE && beta <= VALUE_INFINITE);
            }

            // Sort the PV lines searched so far and update the GUI
            std::stable_sort(rootMoves.begin() + pvFirst, rootMoves.begin() + pvIdx + 1);

            if (mainThread
                && (threads.stop || pvIdx + 1 == multiPV || nodes > 10000000)
                // A thread that aborted search can have mated-in/TB-loss PV and
                // score that cannot be trusted, i.e. it can be delayed or refuted
                // if we would have had time to fully search other root-moves. Thus
                // we suppress this output and below pick a proven score/PV for this
                // thread (from the previous iteration).
                && !(threads.abortedSearch && is_loss(rootMoves[0].uciScore)))
                main_manager()->pv(*this, threads, tt, rootDepth);

            if (threads.stop)
                break;
        }

        if (!threads.stop)
            completedDepth = rootDepth;

        // We make sure not to pick an unproven mated-in score,
        // in case this thread prematurely stopped search (aborted-search).
        if (threads.abortedSearch && rootMoves[0].score != -VALUE_INFINITE
            && is_loss(rootMoves[0].score))
        {
            // Bring the last best move to the front for best thread selection.
            Utility::move_to_front(rootMoves, [&lastBestPV = std::as_const(lastBestPV)](
                                                const auto& rm) { return rm == lastBestPV[0]; });
            rootMoves[0].pv    = lastBestPV;
            rootMoves[0].score = rootMoves[0].uciScore = lastBestScore;
        }
        else if (rootMoves[0].pv[0] != lastBestPV[0])
        {
            lastBestPV        = rootMoves[0].pv;
            lastBestScore     = rootMoves[0].score;
            lastBestMoveDepth = rootDepth;
        }

        if (!mainThread)
            continue;

        // Have we found a "mate in x"?
        if (limits.mate && rootMoves[0].score == rootMoves[0].uciScore
            && ((rootMoves[0].score >= VALUE_MATE_IN_MAX_PLY
                 && VALUE_MATE - rootMoves[0].score <= 2 * limits.mate)
                || (rootMoves[0].score != -VALUE_INFINITE
                    && rootMoves[0].score <= VALUE_MATED_IN_MAX_PLY
                    && VALUE_MATE + rootMoves[0].score <= 2 * limits.mate)))
            threads.stop = true;

        // If the skill level is enabled and time is up, pick a sub-optimal best move
        if (skill.enabled() && skill.time_to_pick(rootDepth))
            skill.pick_best(rootMoves, multiPV);

        // Use part of the gained time from a previous stable move for the current move
        for (auto&& th : threads)
        {
            totBestMoveChanges += th->worker->bestMoveChanges;
            th->worker->bestMoveChanges = 0;
        }

        // Do we have time for the next iteration? Can we stop searching now?
        if (limits.use_time_management() && !threads.stop && !mainThread->stopOnPonderhit)
        {
            int nodesEffort = rootMoves[0].effort * 100000 / std::max(size_t(1), size_t(nodes));

            double fallingEval =
              (11.396 + 2.035 * (mainThread->bestPreviousAverageScore - bestValue)
               + 0.968 * (mainThread->iterValue[iterIdx] - bestValue))
              / 100.0;
            fallingEval = std::clamp(fallingEval, 0.5786, 1.6752);

            // If the bestMove is stable over several iterations, reduce time accordingly
            timeReduction = lastBestMoveDepth + 8 < completedDepth ? 1.4857 : 0.7046;
            double reduction =
              (1.4540 + mainThread->previousTimeReduction) / (2.1593 * timeReduction);
            double bestMoveInstability = 0.9929 + 1.8519 * totBestMoveChanges / threads.size();

            double totalTime =
              mainThread->tm.optimum() * fallingEval * reduction * bestMoveInstability;

            // Cap used time in case of a single legal move for a better viewer experience
            if (rootMoves.size() == 1)
                totalTime = std::min(500.0, totalTime);

            auto elapsedTime = elapsed();

            if (completedDepth >= 10 && nodesEffort >= 97056 && elapsedTime > totalTime * 0.6540
                && !mainThread->ponder)
                threads.stop = true;

            // Stop the search if we have exceeded the totalTime
            if (elapsedTime > totalTime)
            {
                // If we are allowed to ponder do not stop the search now but
                // keep pondering until the GUI sends "ponderhit" or "stop".
                if (mainThread->ponder)
                    mainThread->stopOnPonderhit = true;
                else
                    threads.stop = true;
            }
            else
                threads.increaseDepth = mainThread->ponder || elapsedTime <= totalTime * 0.5138;
        }

        mainThread->iterValue[iterIdx] = bestValue;
        iterIdx                        = (iterIdx + 1) & 3;
    }

    if (!mainThread)
        return;

    mainThread->previousTimeReduction = timeReduction;

    // If the skill level is enabled, swap the best PV line with the sub-optimal one
    if (skill.enabled())
        std::swap(rootMoves[0],
                  *std::find(rootMoves.begin(), rootMoves.end(),
                             skill.best ? skill.best : skill.pick_best(rootMoves, multiPV)));
}

// Reset histories, usually before a new game
void Search::Worker::clear() {
    mainHistory.fill(65);
    lowPlyHistory.fill(107);
    captureHistory.fill(-655);
    pawnHistory.fill(-1215);
    pawnCorrectionHistory.fill(4);
    minorPieceCorrectionHistory.fill(0);
    nonPawnCorrectionHistory[WHITE].fill(0);
    nonPawnCorrectionHistory[BLACK].fill(0);

    for (auto& to : continuationCorrectionHistory)
        for (auto& h : to)
            h.fill(0);

    for (bool inCheck : {false, true})
        for (StatsType c : {NoCaptures, Captures})
            for (auto& to : continuationHistory[inCheck][c])
                for (auto& h : to)
                    h.fill(-493);

    for (size_t i = 1; i < reductions.size(); ++i)
        reductions[i] = int(2937 / 128.0 * std::log(i));

    refreshTable.clear(networks[numaAccessToken]);
}


// Main search function for both PV and non-PV nodes
template<NodeType nodeType>
Value Search::Worker::search(
  Position& pos, Stack* ss, Value alpha, Value beta, Depth depth, bool cutNode) {

    constexpr bool PvNode   = nodeType != NonPV;
    constexpr bool rootNode = nodeType == Root;
    const bool     allNode  = !(PvNode || cutNode);

    // Dive into quiescence search when the depth reaches zero
    if (depth <= 0)
    {
        constexpr auto nt = PvNode ? PV : NonPV;
        return qsearch<nt>(pos, ss, alpha, beta);
    }

    // Limit the depth if extensions made it too large
    depth = std::min(depth, MAX_PLY - 1);

    // Check if we have an upcoming move that draws by repetition
    if (!rootNode && alpha < VALUE_DRAW && pos.upcoming_repetition(ss->ply))
    {
        alpha = value_draw(this->nodes);
        if (alpha >= beta)
            return alpha;
    }

    assert(-VALUE_INFINITE <= alpha && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - 1));
    assert(0 < depth && depth < MAX_PLY);
    assert(!(PvNode && cutNode));

    Move      pv[MAX_PLY + 1];
    StateInfo st;
    ASSERT_ALIGNED(&st, Eval::NNUE::CacheLineSize);

    Key   posKey;
    Move  move, excludedMove, bestMove;
    Depth extension, newDepth;
    Value bestValue, value, eval, maxValue, probCutBeta;
    bool  givesCheck, improving, priorCapture, opponentWorsening;
    bool  capture, ttCapture;
    int   priorReduction = (ss - 1)->reduction;
    (ss - 1)->reduction  = 0;
    Piece movedPiece;

    ValueList<Move, 32> capturesSearched;
    ValueList<Move, 32> quietsSearched;

    // Step 1. Initialize node
    Worker* thisThread = this;
    ss->inCheck        = pos.checkers();
    priorCapture       = pos.captured_piece();
    Color us           = pos.side_to_move();
    ss->moveCount      = 0;
    bestValue          = -VALUE_INFINITE;
    maxValue           = VALUE_INFINITE;

    // Check for the available remaining time
    if (is_mainthread())
        main_manager()->check_time(*thisThread);

    // Used to send selDepth info to GUI (selDepth counts from 1, ply from 0)
    if (PvNode && thisThread->selDepth < ss->ply + 1)
        thisThread->selDepth = ss->ply + 1;

    if (!rootNode)
    {
        // Step 2. Check for aborted search and immediate draw
        if (threads.stop.load(std::memory_order_relaxed) || pos.is_draw(ss->ply)
            || ss->ply >= MAX_PLY)
            return (ss->ply >= MAX_PLY && !ss->inCheck) ? evaluate(pos)
                                                        : value_draw(thisThread->nodes);

        // Step 3. Mate distance pruning. Even if we mate at the next move our score
        // would be at best mate_in(ss->ply + 1), but if alpha is already bigger because
        // a shorter mate was found upward in the tree then there is no need to search
        // because we will never beat the current alpha. Same logic but with reversed
        // signs apply also in the opposite condition of being mated instead of giving
        // mate. In this case, return a fail-high score.
        alpha = std::max(mated_in(ss->ply), alpha);
        beta  = std::min(mate_in(ss->ply + 1), beta);
        if (alpha >= beta)
            return alpha;
    }

    assert(0 <= ss->ply && ss->ply < MAX_PLY);

    bestMove            = Move::none();
    (ss + 2)->cutoffCnt = 0;
    Square prevSq = ((ss - 1)->currentMove).is_ok() ? ((ss - 1)->currentMove).to_sq() : SQ_NONE;
    ss->statScore = 0;

    // Step 4. Transposition table lookup
    excludedMove                   = ss->excludedMove;
    posKey                         = pos.key();
    auto [ttHit, ttData, ttWriter] = tt.probe(posKey);
    // Need further processing of the saved data
    ss->ttHit    = ttHit;
    ttData.move  = rootNode ? thisThread->rootMoves[thisThread->pvIdx].pv[0]
                 : ttHit    ? ttData.move
                            : Move::none();
    ttData.value = ttHit ? value_from_tt(ttData.value, ss->ply, pos.rule50_count()) : VALUE_NONE;
    ss->ttPv     = excludedMove ? ss->ttPv : PvNode || (ttHit && ttData.is_pv);
    ttCapture    = ttData.move && pos.capture_stage(ttData.move);

    // At this point, if excluded, skip straight to step 6, static eval. However,
    // to save indentation, we list the condition in all code between here and there.

    // At non-PV nodes we check for an early TT cutoff
    if (!PvNode && !excludedMove && ttData.depth > depth - (ttData.value <= beta)
        && is_valid(ttData.value)  // Can happen when !ttHit or when access race in probe()
        && (ttData.bound & (ttData.value >= beta ? BOUND_LOWER : BOUND_UPPER))
        && (cutNode == (ttData.value >= beta) || (depth > 9 || (rootDepth > 10 && depth > 5))))
    {
        // If ttMove is quiet, update move sorting heuristics on TT hit
        if (ttData.move && ttData.value >= beta)
        {
            // Bonus for a quiet ttMove that fails high
            if (!ttCapture)
                update_quiet_histories(pos, ss, *this, ttData.move, stat_bonus(depth) * 784 / 1024);

            // Extra penalty for early quiet moves of the previous ply
            if (prevSq != SQ_NONE && (ss - 1)->moveCount <= 3 && !priorCapture)
                update_continuation_histories(ss - 1, pos.piece_on(prevSq), prevSq,
                                              -stat_malus(depth + 1) * 1018 / 1024);
        }

        // Partial workaround for the graph history interaction problem
        // For high rule50 counts don't produce transposition table cutoffs.
        if (pos.rule50_count() < 90)
            return ttData.value;
    }

    // Step 5. Tablebases probe
    if (!rootNode && !excludedMove && tbConfig.cardinality)
    {
        int piecesCount = pos.count<ALL_PIECES>();

        if (piecesCount <= tbConfig.cardinality
            && (piecesCount < tbConfig.cardinality || depth >= tbConfig.probeDepth)
            && pos.rule50_count() == 0 && !pos.can_castle(ANY_CASTLING))
        {
            TB::ProbeState err;
            TB::WDLScore   wdl = Tablebases::probe_wdl(pos, &err);

            // Force check of time on the next occasion
            if (is_mainthread())
                main_manager()->callsCnt = 0;

            if (err != TB::ProbeState::FAIL)
            {
                thisThread->tbHits.fetch_add(1, std::memory_order_relaxed);

                int drawScore = tbConfig.useRule50 ? 1 : 0;

                Value tbValue = VALUE_TB - ss->ply;

                // Use the range VALUE_TB to VALUE_TB_WIN_IN_MAX_PLY to score
                value = wdl < -drawScore ? -tbValue
                      : wdl > drawScore  ? tbValue
                                         : VALUE_DRAW + 2 * wdl * drawScore;

                Bound b = wdl < -drawScore ? BOUND_UPPER
                        : wdl > drawScore  ? BOUND_LOWER
                                           : BOUND_EXACT;

                if (b == BOUND_EXACT || (b == BOUND_LOWER ? value >= beta : value <= alpha))
                {
                    ttWriter.write(posKey, value_to_tt(value, ss->ply), ss->ttPv, b,
                                   std::min(MAX_PLY - 1, depth + 6), Move::none(), VALUE_NONE,
                                   tt.generation());

                    return value;
                }

                if (PvNode)
                {
                    if (b == BOUND_LOWER)
                        bestValue = value, alpha = std::max(alpha, bestValue);
                    else
                        maxValue = value;
                }
            }
        }
    }

    // Step 6. Static evaluation of the position
    Value      unadjustedStaticEval = VALUE_NONE;
    const auto correctionValue      = correction_value(*thisThread, pos, ss);
    if (ss->inCheck)
    {
        // Skip early pruning when in check
        ss->staticEval = eval = (ss - 2)->staticEval;
        improving             = false;
        goto moves_loop;
    }
    else if (excludedMove)
    {
        // Providing the hint that this node's accumulator will be used often
        Eval::NNUE::hint_common_parent_position(pos, networks[numaAccessToken], refreshTable);
        unadjustedStaticEval = eval = ss->staticEval;
    }
    else if (ss->ttHit)
    {
        // Never assume anything about values stored in TT
        unadjustedStaticEval = ttData.eval;
        if (!is_valid(unadjustedStaticEval))
            unadjustedStaticEval = evaluate(pos);
        else if (PvNode)
            Eval::NNUE::hint_common_parent_position(pos, networks[numaAccessToken], refreshTable);

        ss->staticEval = eval = to_corrected_static_eval(unadjustedStaticEval, correctionValue);

        // ttValue can be used as a better position evaluation
        if (is_valid(ttData.value)
            && (ttData.bound & (ttData.value > eval ? BOUND_LOWER : BOUND_UPPER)))
            eval = ttData.value;
    }
    else
    {
        unadjustedStaticEval = evaluate(pos);
        ss->staticEval = eval = to_corrected_static_eval(unadjustedStaticEval, correctionValue);

        // Static evaluation is saved as it was before adjustment by correction history
        ttWriter.write(posKey, VALUE_NONE, ss->ttPv, BOUND_NONE, DEPTH_UNSEARCHED, Move::none(),
                       unadjustedStaticEval, tt.generation());
    }

    // Use static evaluation difference to improve quiet move ordering
    if (((ss - 1)->currentMove).is_ok() && !(ss - 1)->inCheck && !priorCapture)
    {
        int bonus = std::clamp(-10 * int((ss - 1)->staticEval + ss->staticEval), -1906, 1450) + 638;
        thisThread->mainHistory[~us][((ss - 1)->currentMove).from_to()] << bonus * 1136 / 1024;
        if (type_of(pos.piece_on(prevSq)) != PAWN && ((ss - 1)->currentMove).type_of() != PROMOTION)
            thisThread->pawnHistory[pawn_structure_index(pos)][pos.piece_on(prevSq)][prevSq]
              << bonus * 1195 / 1024;
    }

    // Set up the improving flag, which is true if current static evaluation is
    // bigger than the previous static evaluation at our turn (if we were in
    // check at our previous move we go back until we weren't in check) and is
    // false otherwise. The improving flag is used in various pruning heuristics.
    improving = ss->staticEval > (ss - 2)->staticEval;

    opponentWorsening = ss->staticEval > -(ss - 1)->staticEval;

    if (priorReduction >= 3 && !opponentWorsening)
        depth++;
    if (priorReduction >= 1 && depth >= 2 && ss->staticEval + (ss - 1)->staticEval > 200)
        depth--;

    // Step 7. Razoring
    // If eval is really low, skip search entirely and return the qsearch value.
    // For PvNodes, we must have a guard against mates being returned.
    if (!PvNode && eval < alpha - 446 - 303 * depth * depth)
        return qsearch<NonPV>(pos, ss, alpha, beta);

    // Step 8. Futility pruning: child node
    // The depth condition is important for mate finding.
    if (!ss->ttPv && depth < 14
        && eval - futility_margin(depth, cutNode && !ss->ttHit, improving, opponentWorsening)
               - (ss - 1)->statScore / 326 + 37 - std::abs(correctionValue) / 132821
             >= beta
        && eval >= beta && (!ttData.move || ttCapture) && !is_loss(beta) && !is_win(eval))
        return beta + (eval - beta) / 3;

    // Step 9. Null move search with verification search
    if (cutNode && (ss - 1)->currentMove != Move::null() && eval >= beta
        && ss->staticEval >= beta - 21 * depth + 455 - 60 * improving && !excludedMove
        && pos.non_pawn_material(us) && ss->ply >= thisThread->nmpMinPly && !is_loss(beta))
    {
        assert(eval - beta >= 0);

        // Null move dynamic reduction based on depth and eval
        Depth R = std::min(int(eval - beta) / 237, 6) + depth / 3 + 5;

        ss->currentMove                   = Move::null();
        ss->continuationHistory           = &thisThread->continuationHistory[0][0][NO_PIECE][0];
        ss->continuationCorrectionHistory = &thisThread->continuationCorrectionHistory[NO_PIECE][0];

        pos.do_null_move(st, tt);

        Value nullValue = -search<NonPV>(pos, ss + 1, -beta, -beta + 1, depth - R, false);

        pos.undo_null_move();

        // Do not return unproven mate or TB scores
        if (nullValue >= beta && !is_win(nullValue))
        {
            if (thisThread->nmpMinPly || depth < 16)
                return nullValue;

            assert(!thisThread->nmpMinPly);  // Recursive verification is not allowed

            // Do verification search at high depths, with null move pruning disabled
            // until ply exceeds nmpMinPly.
            thisThread->nmpMinPly = ss->ply + 3 * (depth - R) / 4;

            Value v = search<NonPV>(pos, ss, beta - 1, beta, depth - R, false);

            thisThread->nmpMinPly = 0;

            if (v >= beta)
                return nullValue;
        }
    }

    improving |= ss->staticEval >= beta + 97;

    // Step 10. Internal iterative reductions
    // For PV nodes without a ttMove as well as for deep enough cutNodes, we decrease depth.
    // (* Scaler) Especially if they make IIR more aggressive.
    if (((PvNode || cutNode) && depth >= 7 - 3 * PvNode) && !ttData.move)
        depth--;

    // Step 11. ProbCut
    // If we have a good enough capture (or queen promotion) and a reduced search
    // returns a value much above beta, we can (almost) safely prune the previous move.
    probCutBeta = beta + 187 - 55 * improving;
    if (depth >= 3
        && !is_decisive(beta)
        // If value from transposition table is lower than probCutBeta, don't attempt
        // probCut there and in further interactions with transposition table cutoff
        // depth is set to depth - 3 because probCut search has depth set to depth - 4
        // but we also do a move before it. So effective depth is equal to depth - 3.
        && !(is_valid(ttData.value) && ttData.value < probCutBeta))
    {
        assert(probCutBeta < VALUE_INFINITE && probCutBeta > beta);

        MovePicker mp(pos, ttData.move, probCutBeta - ss->staticEval, &thisThread->captureHistory);
        Depth      probCutDepth = std::max(depth - 4, 0);

        while ((move = mp.next_move()) != Move::none())
        {
            assert(move.is_ok());

            if (move == excludedMove)
                continue;

            if (!pos.legal(move))
                continue;

            assert(pos.capture_stage(move));

            movedPiece = pos.moved_piece(move);

            pos.do_move(move, st, &tt);
            thisThread->nodes.fetch_add(1, std::memory_order_relaxed);

            ss->currentMove = move;
            ss->isTTMove    = (move == ttData.move);
            ss->continuationHistory =
              &this->continuationHistory[ss->inCheck][true][movedPiece][move.to_sq()];
            ss->continuationCorrectionHistory =
              &this->continuationCorrectionHistory[movedPiece][move.to_sq()];

            // Perform a preliminary qsearch to verify that the move holds
            value = -qsearch<NonPV>(pos, ss + 1, -probCutBeta, -probCutBeta + 1);

            // If the qsearch held, perform the regular search
            if (value >= probCutBeta && probCutDepth > 0)
                value = -search<NonPV>(pos, ss + 1, -probCutBeta, -probCutBeta + 1, probCutDepth,
                                       !cutNode);

            pos.undo_move(move);

            if (value >= probCutBeta)
            {
                // Save ProbCut data into transposition table
                ttWriter.write(posKey, value_to_tt(value, ss->ply), ss->ttPv, BOUND_LOWER,
                               probCutDepth + 1, move, unadjustedStaticEval, tt.generation());

                if (!is_decisive(value))
                    return value - (probCutBeta - beta);
            }
        }
    }

moves_loop:  // When in check, search starts here

    // Step 12. A small Probcut idea
    probCutBeta = beta + 413;
    if ((ttData.bound & BOUND_LOWER) && ttData.depth >= depth - 4 && ttData.value >= probCutBeta
        && !is_decisive(beta) && is_valid(ttData.value) && !is_decisive(ttData.value))
        return probCutBeta;

    const PieceToHistory* contHist[] = {
      (ss - 1)->continuationHistory, (ss - 2)->continuationHistory, (ss - 3)->continuationHistory,
      (ss - 4)->continuationHistory, (ss - 5)->continuationHistory, (ss - 6)->continuationHistory};


    MovePicker mp(pos, ttData.move, depth, &thisThread->mainHistory, &thisThread->lowPlyHistory,
                  &thisThread->captureHistory, contHist, &thisThread->pawnHistory, ss->ply);

    value = bestValue;

    int moveCount = 0;

    // Step 13. Loop through all pseudo-legal moves until no moves remain
    // or a beta cutoff occurs.
    while ((move = mp.next_move()) != Move::none())
    {
        assert(move.is_ok());

        if (move == excludedMove)
            continue;

        // Check for legality
        if (!pos.legal(move))
            continue;

        // At root obey the "searchmoves" option and skip moves not listed in Root
        // Move List. In MultiPV mode we also skip PV moves that have been already
        // searched and those of lower "TB rank" if we are in a TB root position.
        if (rootNode
            && !std::count(thisThread->rootMoves.begin() + thisThread->pvIdx,
                           thisThread->rootMoves.begin() + thisThread->pvLast, move))
            continue;

        ss->moveCount = ++moveCount;

        if (rootNode && is_mainthread() && nodes > 10000000)
        {
            main_manager()->updates.onIter(
              {depth, UCIEngine::move(move, pos.is_chess960()), moveCount + thisThread->pvIdx});
        }
        if (PvNode)
            (ss + 1)->pv = nullptr;

        extension  = 0;
        capture    = pos.capture_stage(move);
        movedPiece = pos.moved_piece(move);
        givesCheck = pos.gives_check(move);

        // Calculate new depth for this move
        newDepth = depth - 1;

        int delta = beta - alpha;

        Depth r = reduction(improving, depth, moveCount, delta);

        r -= 32 * moveCount;

        // Increase reduction for ttPv nodes (*Scaler)
        // Smaller or even negative value is better for short time controls
        // Bigger value is better for long time controls
        if (ss->ttPv)
            r += 1031;

        // Step 14. Pruning at shallow depth.
        // Depth conditions are important for mate finding.
        if (!rootNode && pos.non_pawn_material(us) && !is_loss(bestValue))
        {
            // Skip quiet moves if movecount exceeds our FutilityMoveCount threshold
            if (moveCount >= futility_move_count(improving, depth))
                mp.skip_quiet_moves();

            // Reduced depth of the next LMR search
            int lmrDepth = newDepth - r / 1024;

            if (capture || givesCheck)
            {
                Piece capturedPiece = pos.piece_on(move.to_sq());
                int   captHist =
                  thisThread->captureHistory[movedPiece][move.to_sq()][type_of(capturedPiece)];

                // Futility pruning for captures
                if (!givesCheck && lmrDepth < 7 && !ss->inCheck)
                {
                    Value futilityValue = ss->staticEval + 242 + 238 * lmrDepth
                                        + PieceValue[capturedPiece] + 95 * captHist / 700;
                    if (futilityValue <= alpha)
                        continue;
                }

                // SEE based pruning for captures and checks
                int seeHist = std::clamp(captHist / 36, -153 * depth, 134 * depth);
                if (!pos.see_ge(move, -157 * depth - seeHist))
                    continue;
            }
            else
            {
                int history =
                  (*contHist[0])[movedPiece][move.to_sq()]
                  + (*contHist[1])[movedPiece][move.to_sq()]
                  + thisThread->pawnHistory[pawn_structure_index(pos)][movedPiece][move.to_sq()];

                // Continuation history based pruning
                if (history < -4107 * depth)
                    continue;

                history += 68 * thisThread->mainHistory[us][move.from_to()] / 32;

                lmrDepth += history / 3576;

                Value futilityValue = ss->staticEval + (bestMove ? 49 : 143) + 116 * lmrDepth;

                if (bestValue < ss->staticEval - 150 && lmrDepth < 7)
                    futilityValue += 108;

                // Futility pruning: parent node
                if (!ss->inCheck && lmrDepth < 12 && futilityValue <= alpha)
                {
                    if (bestValue <= futilityValue && !is_decisive(bestValue)
                        && !is_win(futilityValue))
                        bestValue = futilityValue;
                    continue;
                }

                lmrDepth = std::max(lmrDepth, 0);

                // Prune moves with negative SEE
                if (!pos.see_ge(move, -26 * lmrDepth * lmrDepth))
                    continue;
            }
        }

        // Step 15. Extensions
        // We take care to not overdo to avoid search getting stuck.
        if (ss->ply < thisThread->rootDepth * 2)
        {
            // Singular extension search. If all moves but one
            // fail low on a search of (alpha-s, beta-s), and just one fails high on
            // (alpha, beta), then that move is singular and should be extended. To
            // verify this we do a reduced search on the position excluding the ttMove
            // and if the result is lower than ttValue minus a margin, then we will
            // extend the ttMove. Recursive singular search is avoided.

            // (* Scaler) Generally, higher singularBeta (i.e closer to ttValue)
            // and lower extension margins scale well.

            if (!rootNode && move == ttData.move && !excludedMove
                && depth >= 5 - (thisThread->completedDepth > 32) + ss->ttPv
                && is_valid(ttData.value) && !is_decisive(ttData.value)
                && (ttData.bound & BOUND_LOWER) && ttData.depth >= depth - 3)
            {
                Value singularBeta  = ttData.value - (55 + 81 * (ss->ttPv && !PvNode)) * depth / 58;
                Depth singularDepth = newDepth / 2;

                ss->excludedMove = move;
                value =
                  search<NonPV>(pos, ss, singularBeta - 1, singularBeta, singularDepth, cutNode);
                ss->excludedMove = Move::none();

                if (value < singularBeta)
                {
                    int corrValAdj1  = std::abs(correctionValue) / 265083;
                    int corrValAdj2  = std::abs(correctionValue) / 253680;
                    int doubleMargin = 267 * PvNode - 181 * !ttCapture - corrValAdj1;
                    int tripleMargin =
                      96 + 282 * PvNode - 250 * !ttCapture + 103 * ss->ttPv - corrValAdj2;

                    extension = 1 + (value < singularBeta - doubleMargin)
                              + (value < singularBeta - tripleMargin);

                    depth++;
                }

                // Multi-cut pruning
                // Our ttMove is assumed to fail high based on the bound of the TT entry,
                // and if after excluding the ttMove with a reduced search we fail high
                // over the original beta, we assume this expected cut-node is not
                // singular (multiple moves fail high), and we can prune the whole
                // subtree by returning a softbound.
                else if (value >= beta && !is_decisive(value))
                    return value;

                // Negative extensions
                // If other moves failed high over (ttValue - margin) without the
                // ttMove on a reduced search, but we cannot do multi-cut because
                // (ttValue - margin) is lower than the original beta, we do not know
                // if the ttMove is singular or can do a multi-cut, so we reduce the
                // ttMove in favor of other moves based on some conditions:

                // If the ttMove is assumed to fail high over current beta
                else if (ttData.value >= beta)
                    extension = -3;

                // If we are on a cutNode but the ttMove is not assumed to fail high
                // over current beta
                else if (cutNode)
                    extension = -2;
            }
        }

        // Step 16. Make the move
        pos.do_move(move, st, givesCheck, &tt);
        thisThread->nodes.fetch_add(1, std::memory_order_relaxed);

        // Add extension to new depth
        newDepth += extension;

        // Update the current move (this must be done after singular extension search)
        ss->currentMove = move;
        ss->isTTMove    = (move == ttData.move);
        ss->continuationHistory =
          &thisThread->continuationHistory[ss->inCheck][capture][movedPiece][move.to_sq()];
        ss->continuationCorrectionHistory =
          &thisThread->continuationCorrectionHistory[movedPiece][move.to_sq()];
        uint64_t nodeCount = rootNode ? uint64_t(nodes) : 0;

        // Decrease reduction for PvNodes (*Scaler)
        if (ss->ttPv)
            r -= 2230 + PvNode * 1013 + (ttData.value > alpha) * 925
               + (ttData.depth >= depth) * (971 + cutNode * 1159);

        // These reduction adjustments have no proven non-linear scaling

        r += 316 - moveCount * 32;

        r -= std::abs(correctionValue) / 31568;

        // Increase reduction for cut nodes
        if (cutNode)
            r += 2608 + 1024 * !ttData.move;

        // Increase reduction if ttMove is a capture but the current move is not a capture
        if (ttCapture && !capture)
            r += 1123 + (depth < 8) * 982;

        // Increase reduction if next ply has a lot of fail high
        if ((ss + 1)->cutoffCnt > 3)
            r += 981 + allNode * 833;

        // For first picked move (ttMove) reduce reduction
        else if (move == ttData.move)
            r -= 1982;

        if (capture)
            ss->statScore =
              688 * int(PieceValue[pos.captured_piece()]) / 100
              + thisThread->captureHistory[movedPiece][move.to_sq()][type_of(pos.captured_piece())]
              - 4653;
        else
            ss->statScore = 2 * thisThread->mainHistory[us][move.from_to()]
                          + (*contHist[0])[movedPiece][move.to_sq()]
                          + (*contHist[1])[movedPiece][move.to_sq()] - 3591;

        // Decrease/increase reduction for moves with a good/bad history
        r -= ss->statScore * 1407 / 16384;

        bool              CC = false;
        std::vector<bool> C  = {
          !PvNode,
          !ss->ttPv,
          capture,      
          !capture,     
          givesCheck,    
          !givesCheck, 
          ss->inCheck,
          !ss->inCheck, 
          priorCapture, 
          !priorCapture, 
          improving,   
          !improving,
          cutNode,
          allNode,
          ttCapture,
          !ttCapture,
          ttData.value > alpha,
          ttData.value <= alpha,
          correctionValue > -2322860,
          correctionValue <= -2322860,
          ttData.depth >= depth,
          ttData.depth < depth,
          ss->statScore > -5902,
          ss->statScore <= -5902,
          depth > 5,
          depth <= 5,
          moveCount > 6,
          moveCount <= 6,
          moveCount==1,
          moveCount!=1,
          move==ttData.move,
          move!=ttData.move,
          bool(ttData.move),
          !ttData.move,
          ttHit,
          !ttHit,
          bool(excludedMove),
          !excludedMove,
          (ss-1)->currentMove==Move::null(),
          (ss-1)->currentMove!=Move::null(),
          ttData.bound == BOUND_LOWER,
          ttData.bound == BOUND_EXACT,
          ttData.bound == BOUND_UPPER,
          ttData.bound != BOUND_LOWER,
          ttData.bound != BOUND_EXACT,
          ttData.bound != BOUND_UPPER,
          (ss+1)->cutoffCnt>3,
          (ss+1)->cutoffCnt<=3,
          type_of(movedPiece)==PAWN,
          type_of(movedPiece)==KNIGHT,
          type_of(movedPiece)==BISHOP,
          type_of(movedPiece)==ROOK,
          type_of(movedPiece)==QUEEN,
          type_of(movedPiece)==KING,
          type_of(movedPiece)!=PAWN,
          type_of(movedPiece)!=KNIGHT,
          type_of(movedPiece)!=BISHOP,
          type_of(movedPiece)!=ROOK,
          type_of(movedPiece)!=QUEEN,
          type_of(movedPiece)!=KING,
          !(ss-1)->ttPv,
          (ss-1)->inCheck,
          !(ss-1)->inCheck,
          (ss-1)->ttHit,
          !(ss-1)->ttHit,
          (ss-1)->statScore > -5902,
          (ss-1)->statScore <= -5902,
        };

        // Step 17. Late moves reduction / extension (LMR)
        if (depth >= 2 && moveCount > 1)
        {
            // In general we want to cap the LMR depth search at newDepth, but when
            // reduction is negative, we allow this move a limited search extension
            // beyond the first move depth.
            // To prevent problems when the max value is less than the min value,
            // std::clamp has been replaced by a more robust implementation.

            CC = true;
            /*
            CC = !ss->ttPv;
            std::vector<bool> C = {
                capture,
                !capture,
                givesCheck,
                !givesCheck,
                ss->inCheck,
                !ss->inCheck,
                priorCapture,
                !priorCapture,
                improving,
                !improving,
                cutNode,
                allNode,
            };
            Cond 0: best 32:. expr=101001100101 score=97.7227% support=208228 good=203486
            Cond 1: best 0:. expr=010101100101 score=99.0441% support=5217004 good=5167135
            Cond 2: best 32:. expr=101001100101 score=97.7227% support=208228 good=203486
            Cond 3: best 0:. expr=010101100101 score=99.0441% support=5217004 good=5167135
            Cond 4: best 2:. expr=010110100001 score=99.0271% support=825129 good=817101
            Cond 5: best 0:. expr=010101100101 score=99.0441% support=5217004 good=5167135
            Cond 6: best 0:. expr=010101100101 score=99.0441% support=5217004 good=5167135
            Cond 7: best 54:. expr=101001010101 score=97.3017% support=155470 good=151275
            Cond 8: best 51:. expr=010100101001 score=97.3246% support=1676937 good=1632073
            Cond 9: best 0:. expr=010101100101 score=99.0441% support=5217004 good=5167135
            Cond 10: best 115:. expr=010101100110 score=96.2213% support=1316109 good=1266377
            Cond 11: best 0:. expr=010101100101 score=99.0441% support=5217004 good=5167135

            best 0:. expr=010101100101 score=99.0441% support=5217004 good=5167135
            best 1:. expr=010100100101 score=99.0418% support=6042133 good=5984236
            best 6:. expr=010000100101 score=98.7969% support=6839873 good=6757584
            best 8:. expr=000101100101 score=98.6977% support=6896387 good=6806572
            best 9:. expr=000100100101 score=98.6766% support=7786353 good=7683306
            best 16:. expr=000000100101 score=98.4963% support=8792631 good=8660416
            best 21:. expr=000100100001 score=98.3721% support=9700321 good=9542407
            best 28:. expr=000000100001 score=98.1464% support=10866148 good=10664729
            best 47:. expr=000000100100 score=97.3876% support=10944860 good=10658939
            best 56:. expr=000101000101 score=97.2799% support=16513768 good=16064583
            best 67:. expr=000100000101 score=97.1322% support=18708015 good=18171501
            best 76:. expr=000000000101 score=96.9748% support=20800077 good=20170834
            best 93:. expr=000101000001 score=96.6648% support=24801339 good=23974169
            best 97:. expr=000100000001 score=96.6124% support=26995586 good=26081087
            best 102:. expr=000001000001 score=96.4251% support=27643107 good=26654905
            best 104:. expr=000000000001 score=96.3945% support=29843015 good=28767022
            best 284:. expr=000101000000 score=92.7628% support=42827181 good=39727705
            best 285:. expr=000100000000 score=92.761% support=46090513 good=42754017
            best 300:. expr=000000000000 score=92.2932% support=50580391 good=46682275
             * */
            /*
            bool CC = !ss->ttPv;
            std::vector<bool> C = {
                capture,
                !capture,
                givesCheck,
                !givesCheck,
                ss->inCheck,
                !ss->inCheck,
                priorCapture,
                //!priorCapture,
                improving,
                !improving,
                //cutNode,
                allNode,
                ttCapture,
                !ttCapture,
            };
             * Cond 0: best 105:. expr=101001101101 score=97.9261% support=200253 good=196100
             * Cond 1: best 0:. expr=011010100110 score=100% support=21 good=21
             * Cond 2: best 0:. expr=011010100110 score=100% support=21 good=21
             * Cond 3: best 4:. expr=010110100110 score=99.5699% support=32786 good=32645
             * Cond 4: best 0:. expr=011010100110 score=100% support=21 good=21
             * Cond 5: best 7:. expr=010101101110 score=99.2628% support=171865 good=170598
             * Cond 6: best 0:. expr=011010100110 score=100% support=21 good=21
             * Cond 7: best 66:. expr=010100110110 score=98.442% support=46982 good=46250
             * Cond 8: best 1:. expr=011010101110 score=100% support=21 good=21
             * Cond 9: best 0:. expr=011010100110 score=100% support=21 good=21
             * Cond 10: best 0:. expr=011010100110 score=100% support=21 good=21
             * Cond 11: best 12:. expr=010101101101 score=99.0367% support=5045139 good=4996537
             *
             * best 0:. expr=011010100110 score=100% support=21 good=21
             * best 2:. expr=010010101110 score=99.5702% support=32807 good=32666
             * best 6:. expr=010100101110 score=99.312% support=204651 good=203243
             * best 8:. expr=010100100110 score=99.1496% support=251633 good=249493
             * best 10:. expr=010101101100 score=99.0441% support=5217004 good=5167135
             * best 11:. expr=010100101100 score=99.0418% support=6042133 good=5984236
             * best 24:. expr=010000101100 score=98.7969% support=6839873 good=6757584
             * best 34:. expr=000100101101 score=98.7056% support=7468683 good=7372011
             * best 36:. expr=000100101100 score=98.6766% support=7786353 good=7683306
             * best 47:. expr=000000101101 score=98.5293% support=8443172 good=8318995
             * best 59:. expr=000000101100 score=98.4963% support=8792631 good=8660416
             * best 68:. expr=000100100101 score=98.3895% support=9323850 good=9173694
             * best 71:. expr=000100100100 score=98.3721% support=9700321 good=9542407
             * best 85:. expr=000000100101 score=98.1671% support=10453489 good=10261883
             * best 87:. expr=000000100100 score=98.1464% support=10866148 good=10664729
             * best 167:. expr=000000101000 score=97.3876% support=10944860 good=10658939
             * best 179:. expr=000101001100 score=97.2799% support=16513768 good=16064583
             * best 201:. expr=000100001100 score=97.1322% support=18708015 good=18171501
             * best 224:. expr=000000001100 score=96.9748% support=20800077 good=20170834
             * best 249:. expr=000101000100 score=96.6648% support=24801339 good=23974169
             * best 256:. expr=000100000100 score=96.6124% support=26995586 good=26081087
             * best 268:. expr=000001000100 score=96.4251% support=27643107 good=26654905
             * best 270:. expr=000000000100 score=96.3945% support=29843015 good=28767022
             * best 557:. expr=000101000001 score=92.7655% support=40574924 good=37639546
             * best 558:. expr=000100000001 score=92.7629% support=43627283 good=40469941
             * best 560:. expr=000100000000 score=92.761% support=46090513 good=42754017
             * best 582:. expr=000000000001 score=92.3158% support=47912101 good=44230425
             * best 585:. expr=000000000000 score=92.2932% support=50580391 good=46682275
             */
            /*
            bool CC = !ss->ttPv;
            std::vector<bool> C = {
                capture,
                !capture,
                givesCheck,
                !givesCheck,
                ss->inCheck,
                !ss->inCheck,
                priorCapture,
                //!priorCapture,
                //improving,
                !improving,
                //cutNode,
                allNode,
                ttCapture,
                //!ttCapture,
                ttData.value > alpha,
                ttData.value <= alpha,
            };
             * Cond 0: best 8:. expr=101010101101 score=100% support=17 good=17
             * Cond 1: best 4:. expr=011010011101 score=100% support=24 good=24
             * Cond 2: best 0:. expr=001010011101 score=100% support=41 good=41
             * Cond 3: best 18:. expr=010110101101 score=99.6351% support=20007 good=19934
             * Cond 4: best 0:. expr=001010011101 score=100% support=41 good=41
             * Cond 5: best 29:. expr=010101111101 score=99.3643% support=105390 good=104720
             * Cond 6: best 2:. expr=001010111101 score=100% support=31 good=31
             * Cond 7: best 0:. expr=001010011101 score=100% support=41 good=41
             * Cond 8: best 0:. expr=001010011101 score=100% support=41 good=41
             * Cond 9: best 0:. expr=001010011101 score=100% support=41 good=41
             * Cond 10: best 14:. expr=011010101110 score=100% support=7 good=7
             * Cond 11: best 0:. expr=001010011101 score=100% support=41 good=41
             *
             * best 0:. expr=001010011101 score=100% support=41 good=41
             * best 16:. expr=010010101101 score=99.6354% support=20021 good=19948
             * best 20:. expr=010010111100 score=99.5702% support=32807 good=32666
             * best 28:. expr=010100111101 score=99.4075% support=125397 good=124654
             * best 34:. expr=010100111100 score=99.312% support=204651 good=203243
             * best 43:. expr=010100101100 score=99.1496% support=251633 good=249493
             * best 48:. expr=010110111001 score=99.081% support=654736 good=648719
             * best 50:. expr=010010101001 score=99.0765% support=655761 good=649705
             * best 52:. expr=010100111001 score=99.0655% support=4573949 good=4531204
             * best 55:. expr=010101111000 score=99.0441% support=5217004 good=5167135
             * best 56:. expr=010100111000 score=99.0418% support=6042133 good=5984236
             * best 90:. expr=010000111000 score=98.7969% support=6839873 good=6757584
             * best 105:. expr=000101111000 score=98.6977% support=6896387 good=6806572
             * best 109:. expr=000100111000 score=98.6766% support=7786353 good=7683306
             * best 144:. expr=000000111000 score=98.4963% support=8792631 good=8660416
             * best 161:. expr=000100101000 score=98.3721% support=9700321 good=9542407
             * best 194:. expr=000000101000 score=98.1464% support=10866148 good=10664729
             * best 327:. expr=000000110000 score=97.3876% support=10944860 good=10658939
             * best 329:. expr=000101011001 score=97.3847% support=11168813 good=10876713
             * best 344:. expr=000101011000 score=97.2799% support=16513768 good=16064583
             * best 369:. expr=000100011000 score=97.1322% support=18708015 good=18171501
             * best 409:. expr=000000011000 score=96.9748% support=20800077 good=20170834
             * best 458:. expr=000101001000 score=96.6648% support=24801339 good=23974169
             * best 465:. expr=000100001000 score=96.6124% support=26995586 good=26081087
             * best 490:. expr=000001001000 score=96.4251% support=27643107 good=26654905
             * best 495:. expr=000000001000 score=96.3945% support=29843015 good=28767022
             * best 898:. expr=000101000000 score=92.7628% support=42827181 good=39727705
             * best 899:. expr=000100000000 score=92.761% support=46090513 good=42754017
             * best 932:. expr=000000000000 score=92.2932% support=50580391 good=46682275
             */
            /*
            bool CC = !ss->ttPv;
            std::vector<bool> C = {
                capture,
                !capture,
                givesCheck,
                //!givesCheck,
                ss->inCheck,
                //!ss->inCheck,
                priorCapture,
                //!priorCapture,
                //improving,
                !improving,
                //cutNode,
                allNode,
                ttCapture,
                //!ttCapture,
                ttData.value > alpha,
                ttData.value <= alpha,
                correctionValue > -2322860,
                correctionValue <= -2322860,
                }:
             * Cond 0: best 14:. expr=101100110100 score=100% support=17 good=17
             * Cond 1: best 4:. expr=011101110100 score=100% support=24 good=24
             * Cond 2: best 0:. expr=001101110100 score=100% support=41 good=41
             * Cond 3: best 0:. expr=001101110100 score=100% support=41 good=41
             * Cond 4: best 2:. expr=001110110100 score=100% support=31 good=31
             * Cond 5: best 0:. expr=001101110100 score=100% support=41 good=41
             * Cond 6: best 0:. expr=001101110100 score=100% support=41 good=41
             * Cond 7: best 0:. expr=001101110100 score=100% support=41 good=41
             * Cond 8: best 43:. expr=011100111010 score=100% support=8 good=8
             * Cond 9: best 0:. expr=001101110100 score=100% support=41 good=41
             * Cond 10: best 10:. expr=001101110110 score=100% support=20 good=20
             * Cond 11: best 7:. expr=001101110101 score=100% support=21 good=21
             *
             * best 0:. expr=001101110100 score=100% support=41 good=41
             * best 54:. expr=010111110101 score=99.7015% support=9044 good=9017
             * best 56:. expr=010110110100 score=99.6354% support=20021 good=19948
             * best 64:. expr=010110110000 score=99.5702% support=32807 good=32666
             * best 78:. expr=010011110101 score=99.2943% support=62491 good=62050
             * best 85:. expr=010010110101 score=99.1855% support=72925 good=72331
             * best 86:. expr=010110100101 score=99.1348% support=280861 good=278431
             * best 90:. expr=010110100100 score=99.0765% support=655761 good=649705
             * best 98:. expr=010110100000 score=99.0229% support=826408 good=818333
             * best 107:. expr=010011100101 score=98.9408% support=2467661 good=2441523
             * best 112:. expr=010011100001 score=98.9128% support=3351317 good=3314880
             * best 126:. expr=010011100100 score=98.8198% support=5199589 good=5138226
             * best 130:. expr=010011100000 score=98.7969% support=6839873 good=6757584
             * best 190:. expr=000011100000 score=98.4963% support=8792631 good=8660416
             * best 240:. expr=000010100000 score=98.1464% support=10866148 good=10664729
             * best 375:. expr=000011000000 score=97.3876% support=10944860 good=10658939
             * best 412:. expr=000001100100 score=97.1396% support=14149864 good=13745128
             * best 442:. expr=000001100000 score=96.9748% support=20800077 good=20170834
             * best 558:. expr=000000100000 score=96.3945% support=29843015 good=28767022
             * best 1112:. expr=000000000000 score=92.2932% support=50580391 good=46682275
             */

            Depth d = std::max(
              1, std::min(newDepth - r / 1024, newDepth + !allNode + (PvNode && !bestMove)));

            ss->reduction = newDepth - d;

            value         = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, d, true);
            ss->reduction = 0;

            if(CC)
            {
                /*
                dbg_mean_of(correctionValue, 0);
                dbg_hit_on(correctionValue > -2322860, 0);
                dbg_mean_of(ss->statScore, 1);
                dbg_hit_on(ss->statScore > -5902, 1);
                dbg_mean_of(depth, 2);
                dbg_hit_on(depth > 5, 2);
                dbg_mean_of(moveCount, 3);
                dbg_hit_on(moveCount > 6, 3);
                */
                //dbg_hit_on(C[0], 0);
                //dbg_hit_on(C[1], 1);
                bool T = value <= alpha;
                learn(T, C);
            }

            // Do a full-depth search when reduced LMR search fails high
            if (value > alpha && d < newDepth)
            {
                // Adjust full-depth search based on LMR results - if the result was
                // good enough search deeper, if it was bad enough search shallower.
                const bool doDeeperSearch    = value > (bestValue + 41 + 2 * newDepth);
                const bool doShallowerSearch = value < bestValue + 9;

                newDepth += doDeeperSearch - doShallowerSearch;

                if (newDepth > d)
                    value = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode);

                // Post LMR continuation history updates
                int bonus = (value >= beta) * 2010;
                update_continuation_histories(ss, movedPiece, move.to_sq(), bonus);
            }
        }

        // Step 18. Full-depth search when LMR is skipped
        else if (!PvNode || moveCount > 1)
        {
            CC = true;

            // Increase reduction if ttMove is not present
            if (!ttData.move)
                r += 1111;

            // Note that if expected reduction is high, we reduce search depth here
            value = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha,
                                   newDepth - (r > 3554) - (r > 5373 && newDepth > 2), !cutNode);

            if (CC)
            {
                /*
                dbg_mean_of(correctionValue, 0);
                dbg_hit_on(correctionValue > -2322860, 0);
                dbg_mean_of(ss->statScore, 1);
                dbg_hit_on(ss->statScore > -5902, 1);
                dbg_mean_of(depth, 2);
                dbg_hit_on(depth > 5, 2);
                dbg_mean_of(moveCount, 3);
                dbg_hit_on(moveCount > 6, 3);
                */
                //dbg_hit_on(C[0], 0);
                //dbg_hit_on(C[1], 1);
                bool T = value <= alpha;
                learn(T, C);
            }
        }

        // For PV nodes only, do a full PV search on the first move or after a fail high,
        // otherwise let the parent node fail low with value <= alpha and try another move.
        if (PvNode && (moveCount == 1 || value > alpha))
        {
            (ss + 1)->pv    = pv;
            (ss + 1)->pv[0] = Move::none();

            // Extend move from transposition table if we are about to dive into qsearch.
            if (move == ttData.move && thisThread->rootDepth > 8)
                newDepth = std::max(newDepth, 1);

            value = -search<PV>(pos, ss + 1, -beta, -alpha, newDepth, false);
        }

        // Step 19. Undo move
        pos.undo_move(move);

        assert(value > -VALUE_INFINITE && value < VALUE_INFINITE);

        // Step 20. Check for a new best move
        // Finished searching the move. If a stop occurred, the return value of
        // the search cannot be trusted, and we return immediately without updating
        // best move, principal variation nor transposition table.
        if (threads.stop.load(std::memory_order_relaxed))
            return VALUE_ZERO;

        if (rootNode)
        {
            RootMove& rm =
              *std::find(thisThread->rootMoves.begin(), thisThread->rootMoves.end(), move);

            rm.effort += nodes - nodeCount;

            rm.averageScore =
              rm.averageScore != -VALUE_INFINITE ? (value + rm.averageScore) / 2 : value;

            rm.meanSquaredScore = rm.meanSquaredScore != -VALUE_INFINITE * VALUE_INFINITE
                                  ? (value * std::abs(value) + rm.meanSquaredScore) / 2
                                  : value * std::abs(value);

            // PV move or new best move?
            if (moveCount == 1 || value > alpha)
            {
                rm.score = rm.uciScore = value;
                rm.selDepth            = thisThread->selDepth;
                rm.scoreLowerbound = rm.scoreUpperbound = false;

                if (value >= beta)
                {
                    rm.scoreLowerbound = true;
                    rm.uciScore        = beta;
                }
                else if (value <= alpha)
                {
                    rm.scoreUpperbound = true;
                    rm.uciScore        = alpha;
                }

                rm.pv.resize(1);

                assert((ss + 1)->pv);

                for (Move* m = (ss + 1)->pv; *m != Move::none(); ++m)
                    rm.pv.push_back(*m);

                // We record how often the best move has been changed in each iteration.
                // This information is used for time management. In MultiPV mode,
                // we must take care to only do this for the first PV line.
                if (moveCount > 1 && !thisThread->pvIdx)
                    ++thisThread->bestMoveChanges;
            }
            else
                // All other moves but the PV, are set to the lowest value: this
                // is not a problem when sorting because the sort is stable and the
                // move position in the list is preserved - just the PV is pushed up.
                rm.score = -VALUE_INFINITE;
        }

        // In case we have an alternative move equal in eval to the current bestmove,
        // promote it to bestmove by pretending it just exceeds alpha (but not beta).
        int inc = (value == bestValue && ss->ply + 2 >= thisThread->rootDepth
                   && (int(nodes) & 15) == 0 && !is_win(std::abs(value) + 1));

        if (value + inc > bestValue)
        {
            bestValue = value;

            if (value + inc > alpha)
            {
                bestMove = move;

                if (PvNode && !rootNode)  // Update pv even in fail-high case
                    update_pv(ss->pv, move, (ss + 1)->pv);

                if (value >= beta)
                {
                    // (* Scaler) Especially if they make cutoffCnt increment more often.
                    ss->cutoffCnt += (extension < 2) || PvNode;
                    assert(value >= beta);  // Fail high
                    break;
                }
                else
                {
                    // Reduce other moves if we have found at least one score improvement
                    if (depth > 2 && depth < 15 && !is_decisive(value))
                        depth -= 2;

                    assert(depth > 0);
                    alpha = value;  // Update alpha! Always alpha < beta
                }
            }
        }

        // If the move is worse than some previously searched move,
        // remember it, to update its stats later.
        if (move != bestMove && moveCount <= 32)
        {
            if (capture)
                capturesSearched.push_back(move);
            else
                quietsSearched.push_back(move);
        }
    }

    // Step 21. Check for mate and stalemate
    // All legal moves have been searched and if there are no legal moves, it
    // must be a mate or a stalemate. If we are in a singular extension search then
    // return a fail low score.

    assert(moveCount || !ss->inCheck || excludedMove || !MoveList<LEGAL>(pos).size());

    // Adjust best value for fail high cases at non-pv nodes
    if (!PvNode && bestValue >= beta && !is_decisive(bestValue) && !is_decisive(beta)
        && !is_decisive(alpha))
        bestValue = (bestValue * depth + beta) / (depth + 1);

    if (!moveCount)
        bestValue = excludedMove ? alpha : ss->inCheck ? mated_in(ss->ply) : VALUE_DRAW;

    // If there is a move that produces search value greater than alpha,
    // we update the stats of searched moves.
    else if (bestMove)
        update_all_stats(pos, ss, *this, bestMove, prevSq, quietsSearched, capturesSearched, depth,
                         bestMove == ttData.move, moveCount);

    // Bonus for prior countermove that caused the fail low
    else if (!priorCapture && prevSq != SQ_NONE)
    {
        int bonusScale = (118 * (depth > 5) + 36 * !allNode + 161 * ((ss - 1)->moveCount > 8)
                          + 133 * (!ss->inCheck && bestValue <= ss->staticEval - 107)
                          + 120 * (!(ss - 1)->inCheck && bestValue <= -(ss - 1)->staticEval - 84)
                          + 81 * ((ss - 1)->isTTMove) + 100 * (ss->cutoffCnt <= 3)
                          + std::min(-(ss - 1)->statScore / 108, 320));

        bonusScale = std::max(bonusScale, 0);

        const int scaledBonus = stat_bonus(depth) * bonusScale;

        update_continuation_histories(ss - 1, pos.piece_on(prevSq), prevSq,
                                      scaledBonus * 416 / 32768);

        thisThread->mainHistory[~us][((ss - 1)->currentMove).from_to()]
          << scaledBonus * 219 / 32768;

        if (type_of(pos.piece_on(prevSq)) != PAWN && ((ss - 1)->currentMove).type_of() != PROMOTION)
            thisThread->pawnHistory[pawn_structure_index(pos)][pos.piece_on(prevSq)][prevSq]
              << scaledBonus * 1103 / 32768;
    }

    else if (priorCapture && prevSq != SQ_NONE)
    {
        // bonus for prior countermoves that caused the fail low
        Piece capturedPiece = pos.captured_piece();
        assert(capturedPiece != NO_PIECE);
        thisThread->captureHistory[pos.piece_on(prevSq)][prevSq][type_of(capturedPiece)]
          << stat_bonus(depth) * 2;
    }

    if (PvNode)
        bestValue = std::min(bestValue, maxValue);

    // If no good move is found and the previous position was ttPv, then the previous
    // opponent move is probably good and the new position is added to the search tree.
    if (bestValue <= alpha)
        ss->ttPv = ss->ttPv || (ss - 1)->ttPv;

    // Write gathered information in transposition table. Note that the
    // static evaluation is saved as it was before correction history.
    if (!excludedMove && !(rootNode && thisThread->pvIdx))
        ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), ss->ttPv,
                       bestValue >= beta    ? BOUND_LOWER
                       : PvNode && bestMove ? BOUND_EXACT
                                            : BOUND_UPPER,
                       depth, bestMove, unadjustedStaticEval, tt.generation());

    // Adjust correction history
    if (!ss->inCheck && !(bestMove && pos.capture(bestMove))
        && ((bestValue < ss->staticEval && bestValue < beta)  // negative correction & no fail high
            || (bestValue > ss->staticEval && bestMove)))     // positive correction & no fail low
    {
        auto bonus = std::clamp(int(bestValue - ss->staticEval) * depth / 8,
                                -CORRECTION_HISTORY_LIMIT / 4, CORRECTION_HISTORY_LIMIT / 4);
        update_correction_history(pos, ss, *thisThread, bonus);
    }

    assert(bestValue > -VALUE_INFINITE && bestValue < VALUE_INFINITE);

    return bestValue;
}


// Quiescence search function, which is called by the main search function with
// depth zero, or recursively with further decreasing depth. With depth <= 0, we
// "should" be using static eval only, but tactical moves may confuse the static eval.
// To fight this horizon effect, we implement this qsearch of tactical moves.
// See https://www.chessprogramming.org/Horizon_Effect
// and https://www.chessprogramming.org/Quiescence_Search
template<NodeType nodeType>
Value Search::Worker::qsearch(Position& pos, Stack* ss, Value alpha, Value beta) {

    static_assert(nodeType != Root);
    constexpr bool PvNode = nodeType == PV;

    assert(alpha >= -VALUE_INFINITE && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - 1));

    // Check if we have an upcoming move that draws by repetition
    if (alpha < VALUE_DRAW && pos.upcoming_repetition(ss->ply))
    {
        alpha = value_draw(this->nodes);
        if (alpha >= beta)
            return alpha;
    }

    Move      pv[MAX_PLY + 1];
    StateInfo st;
    ASSERT_ALIGNED(&st, Eval::NNUE::CacheLineSize);

    Key   posKey;
    Move  move, bestMove;
    Value bestValue, value, futilityBase;
    bool  pvHit, givesCheck, capture;
    int   moveCount;

    // Step 1. Initialize node
    if (PvNode)
    {
        (ss + 1)->pv = pv;
        ss->pv[0]    = Move::none();
    }

    Worker* thisThread = this;
    bestMove           = Move::none();
    ss->inCheck        = pos.checkers();
    moveCount          = 0;

    // Used to send selDepth info to GUI (selDepth counts from 1, ply from 0)
    if (PvNode && thisThread->selDepth < ss->ply + 1)
        thisThread->selDepth = ss->ply + 1;

    // Step 2. Check for an immediate draw or maximum ply reached
    if (pos.is_draw(ss->ply) || ss->ply >= MAX_PLY)
        return (ss->ply >= MAX_PLY && !ss->inCheck) ? evaluate(pos) : VALUE_DRAW;

    assert(0 <= ss->ply && ss->ply < MAX_PLY);

    // Step 3. Transposition table lookup
    posKey                         = pos.key();
    auto [ttHit, ttData, ttWriter] = tt.probe(posKey);
    // Need further processing of the saved data
    ss->ttHit    = ttHit;
    ttData.move  = ttHit ? ttData.move : Move::none();
    ttData.value = ttHit ? value_from_tt(ttData.value, ss->ply, pos.rule50_count()) : VALUE_NONE;
    pvHit        = ttHit && ttData.is_pv;

    // At non-PV nodes we check for an early TT cutoff
    if (!PvNode && ttData.depth >= DEPTH_QS
        && is_valid(ttData.value)  // Can happen when !ttHit or when access race in probe()
        && (ttData.bound & (ttData.value >= beta ? BOUND_LOWER : BOUND_UPPER)))
        return ttData.value;

    // Step 4. Static evaluation of the position
    Value      unadjustedStaticEval = VALUE_NONE;
    const auto correctionValue      = correction_value(*thisThread, pos, ss);
    if (ss->inCheck)
        bestValue = futilityBase = -VALUE_INFINITE;
    else
    {
        if (ss->ttHit)
        {
            // Never assume anything about values stored in TT
            unadjustedStaticEval = ttData.eval;
            if (!is_valid(unadjustedStaticEval))
                unadjustedStaticEval = evaluate(pos);
            ss->staticEval = bestValue =
              to_corrected_static_eval(unadjustedStaticEval, correctionValue);

            // ttValue can be used as a better position evaluation
            if (is_valid(ttData.value) && !is_decisive(ttData.value)
                && (ttData.bound & (ttData.value > bestValue ? BOUND_LOWER : BOUND_UPPER)))
                bestValue = ttData.value;
        }
        else
        {
            // In case of null move search, use previous static eval with opposite sign
            unadjustedStaticEval =
              (ss - 1)->currentMove != Move::null() ? evaluate(pos) : -(ss - 1)->staticEval;
            ss->staticEval = bestValue =
              to_corrected_static_eval(unadjustedStaticEval, correctionValue);
        }

        // Stand pat. Return immediately if static value is at least beta
        if (bestValue >= beta)
        {
            if (!is_decisive(bestValue))
                bestValue = (bestValue + beta) / 2;
            if (!ss->ttHit)
                ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), false, BOUND_LOWER,
                               DEPTH_UNSEARCHED, Move::none(), unadjustedStaticEval,
                               tt.generation());
            return bestValue;
        }

        if (bestValue > alpha)
            alpha = bestValue;

        futilityBase = ss->staticEval + 325;
    }

    const PieceToHistory* contHist[] = {(ss - 1)->continuationHistory,
                                        (ss - 2)->continuationHistory};

    Square prevSq = ((ss - 1)->currentMove).is_ok() ? ((ss - 1)->currentMove).to_sq() : SQ_NONE;

    // Initialize a MovePicker object for the current position, and prepare to search
    // the moves. We presently use two stages of move generator in quiescence search:
    // captures, or evasions only when in check.
    MovePicker mp(pos, ttData.move, DEPTH_QS, &thisThread->mainHistory, &thisThread->lowPlyHistory,
                  &thisThread->captureHistory, contHist, &thisThread->pawnHistory, ss->ply);

    // Step 5. Loop through all pseudo-legal moves until no moves remain or a beta
    // cutoff occurs.
    while ((move = mp.next_move()) != Move::none())
    {
        assert(move.is_ok());

        if (!pos.legal(move))
            continue;

        givesCheck = pos.gives_check(move);
        capture    = pos.capture_stage(move);

        moveCount++;

        // Step 6. Pruning
        if (!is_loss(bestValue))
        {
            // Futility pruning and moveCount pruning
            if (!givesCheck && move.to_sq() != prevSq && !is_loss(futilityBase)
                && move.type_of() != PROMOTION)
            {
                if (moveCount > 2)
                    continue;

                Value futilityValue = futilityBase + PieceValue[pos.piece_on(move.to_sq())];

                // If static eval + value of piece we are going to capture is
                // much lower than alpha, we can prune this move.
                if (futilityValue <= alpha)
                {
                    bestValue = std::max(bestValue, futilityValue);
                    continue;
                }

                // If static exchange evaluation is low enough
                // we can prune this move.
                if (!pos.see_ge(move, alpha - futilityBase))
                {
                    bestValue = std::min(alpha, futilityBase);
                    continue;
                }
            }

            // Continuation history based pruning
            if (!capture
                && (*contHist[0])[pos.moved_piece(move)][move.to_sq()]
                       + (*contHist[1])[pos.moved_piece(move)][move.to_sq()]
                       + thisThread->pawnHistory[pawn_structure_index(pos)][pos.moved_piece(move)]
                                                [move.to_sq()]
                     <= 5389)
                continue;

            // Do not search moves with bad enough SEE values
            if (!pos.see_ge(move, -75))
                continue;
        }

        // Step 7. Make and search the move
        Piece movedPiece = pos.moved_piece(move);

        pos.do_move(move, st, givesCheck, &tt);
        thisThread->nodes.fetch_add(1, std::memory_order_relaxed);

        // Update the current move
        ss->currentMove = move;
        ss->continuationHistory =
          &thisThread->continuationHistory[ss->inCheck][capture][movedPiece][move.to_sq()];
        ss->continuationCorrectionHistory =
          &thisThread->continuationCorrectionHistory[movedPiece][move.to_sq()];

        value = -qsearch<nodeType>(pos, ss + 1, -beta, -alpha);
        pos.undo_move(move);

        assert(value > -VALUE_INFINITE && value < VALUE_INFINITE);

        // Step 8. Check for a new best move
        if (value > bestValue)
        {
            bestValue = value;

            if (value > alpha)
            {
                bestMove = move;

                if (PvNode)  // Update pv even in fail-high case
                    update_pv(ss->pv, move, (ss + 1)->pv);

                if (value < beta)  // Update alpha here!
                    alpha = value;
                else
                    break;  // Fail high
            }
        }
    }

    // Step 9. Check for mate
    // All legal moves have been searched. A special case: if we are
    // in check and no legal moves were found, it is checkmate.
    if (ss->inCheck && bestValue == -VALUE_INFINITE)
    {
        assert(!MoveList<LEGAL>(pos).size());
        return mated_in(ss->ply);  // Plies to mate from the root
    }

    if (!is_decisive(bestValue) && bestValue >= beta)
        bestValue = (3 * bestValue + beta) / 4;

    // Save gathered info in transposition table. The static evaluation
    // is saved as it was before adjustment by correction history.
    ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), pvHit,
                   bestValue >= beta ? BOUND_LOWER : BOUND_UPPER, DEPTH_QS, bestMove,
                   unadjustedStaticEval, tt.generation());

    assert(bestValue > -VALUE_INFINITE && bestValue < VALUE_INFINITE);

    return bestValue;
}

Depth Search::Worker::reduction(bool i, Depth d, int mn, int delta) const {
    int reductionScale = reductions[d] * reductions[mn];
    return reductionScale - delta * 735 / rootDelta + !i * reductionScale * 191 / 512 + 1132;
}

// elapsed() returns the time elapsed since the search started. If the
// 'nodestime' option is enabled, it will return the count of nodes searched
// instead. This function is called to check whether the search should be
// stopped based on predefined thresholds like time limits or nodes searched.
//
// elapsed_time() returns the actual time elapsed since the start of the search.
// This function is intended for use only when printing PV outputs, and not used
// for making decisions within the search algorithm itself.
TimePoint Search::Worker::elapsed() const {
    return main_manager()->tm.elapsed([this]() { return threads.nodes_searched(); });
}

TimePoint Search::Worker::elapsed_time() const { return main_manager()->tm.elapsed_time(); }

Value Search::Worker::evaluate(const Position& pos) {
    return Eval::evaluate(networks[numaAccessToken], pos, refreshTable,
                          optimism[pos.side_to_move()]);
}

namespace {
// Adjusts a mate or TB score from "plies to mate from the root" to
// "plies to mate from the current position". Standard scores are unchanged.
// The function is called before storing a value in the transposition table.
Value value_to_tt(Value v, int ply) { return is_win(v) ? v + ply : is_loss(v) ? v - ply : v; }


// Inverse of value_to_tt(): it adjusts a mate or TB score from the transposition
// table (which refers to the plies to mate/be mated from current position) to
// "plies to mate/be mated (TB win/loss) from the root". However, to avoid
// potentially false mate or TB scores related to the 50 moves rule and the
// graph history interaction, we return the highest non-TB score instead.
Value value_from_tt(Value v, int ply, int r50c) {

    if (!is_valid(v))
        return VALUE_NONE;

    // handle TB win or better
    if (is_win(v))
    {
        // Downgrade a potentially false mate score
        if (v >= VALUE_MATE_IN_MAX_PLY && VALUE_MATE - v > 100 - r50c)
            return VALUE_TB_WIN_IN_MAX_PLY - 1;

        // Downgrade a potentially false TB score.
        if (VALUE_TB - v > 100 - r50c)
            return VALUE_TB_WIN_IN_MAX_PLY - 1;

        return v - ply;
    }

    // handle TB loss or worse
    if (is_loss(v))
    {
        // Downgrade a potentially false mate score.
        if (v <= VALUE_MATED_IN_MAX_PLY && VALUE_MATE + v > 100 - r50c)
            return VALUE_TB_LOSS_IN_MAX_PLY + 1;

        // Downgrade a potentially false TB score.
        if (VALUE_TB + v > 100 - r50c)
            return VALUE_TB_LOSS_IN_MAX_PLY + 1;

        return v + ply;
    }

    return v;
}


// Adds current move and appends child pv[]
void update_pv(Move* pv, Move move, const Move* childPv) {

    for (*pv++ = move; childPv && *childPv != Move::none();)
        *pv++ = *childPv++;
    *pv = Move::none();
}


// Updates stats at the end of search() when a bestMove is found
void update_all_stats(const Position&      pos,
                      Stack*               ss,
                      Search::Worker&      workerThread,
                      Move                 bestMove,
                      Square               prevSq,
                      ValueList<Move, 32>& quietsSearched,
                      ValueList<Move, 32>& capturesSearched,
                      Depth                depth,
                      bool                 isTTMove,
                      int                  moveCount) {

    CapturePieceToHistory& captureHistory = workerThread.captureHistory;
    Piece                  moved_piece    = pos.moved_piece(bestMove);
    PieceType              captured;

    int bonus = stat_bonus(depth) + 298 * isTTMove;
    int malus = stat_malus(depth) - 32 * (moveCount - 1);

    if (!pos.capture_stage(bestMove))
    {
        update_quiet_histories(pos, ss, workerThread, bestMove, bonus * 1202 / 1024);

        // Decrease stats for all non-best quiet moves
        for (Move move : quietsSearched)
            update_quiet_histories(pos, ss, workerThread, move, -malus * 1152 / 1024);
    }
    else
    {
        // Increase stats for the best move in case it was a capture move
        captured = type_of(pos.piece_on(bestMove.to_sq()));
        captureHistory[moved_piece][bestMove.to_sq()][captured] << bonus * 1236 / 1024;
    }

    // Extra penalty for a quiet early move that was not a TT move in
    // previous ply when it gets refuted.
    if (prevSq != SQ_NONE && ((ss - 1)->moveCount == 1 + (ss - 1)->ttHit) && !pos.captured_piece())
        update_continuation_histories(ss - 1, pos.piece_on(prevSq), prevSq, -malus * 976 / 1024);

    // Decrease stats for all non-best capture moves
    for (Move move : capturesSearched)
    {
        moved_piece = pos.moved_piece(move);
        captured    = type_of(pos.piece_on(move.to_sq()));
        captureHistory[moved_piece][move.to_sq()][captured] << -malus * 1224 / 1024;
    }
}


// Updates histories of the move pairs formed by moves
// at ply -1, -2, -3, -4, and -6 with current move.
void update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus) {
    static constexpr std::array<ConthistBonus, 6> conthist_bonuses = {
      {{1, 1029}, {2, 656}, {3, 326}, {4, 536}, {5, 120}, {6, 537}}};

    for (const auto [i, weight] : conthist_bonuses)
    {
        // Only update the first 2 continuation histories if we are in check
        if (ss->inCheck && i > 2)
            break;
        if (((ss - i)->currentMove).is_ok())
            (*(ss - i)->continuationHistory)[pc][to] << bonus * weight / 1024;
    }
}

// Updates move sorting heuristics

void update_quiet_histories(
  const Position& pos, Stack* ss, Search::Worker& workerThread, Move move, int bonus) {

    Color us = pos.side_to_move();
    workerThread.mainHistory[us][move.from_to()] << bonus;  // Untuned to prevent duplicate effort

    if (ss->ply < LOW_PLY_HISTORY_SIZE)
        workerThread.lowPlyHistory[ss->ply][move.from_to()] << bonus * 844 / 1024;

    update_continuation_histories(ss, pos.moved_piece(move), move.to_sq(), bonus * 964 / 1024);

    int pIndex = pawn_structure_index(pos);
    workerThread.pawnHistory[pIndex][pos.moved_piece(move)][move.to_sq()] << bonus * 615 / 1024;
}

}

// When playing with strength handicap, choose the best move among a set of
// RootMoves using a statistical rule dependent on 'level'. Idea by Heinz van Saanen.
Move Skill::pick_best(const RootMoves& rootMoves, size_t multiPV) {
    static PRNG rng(now());  // PRNG sequence should be non-deterministic

    // RootMoves are already sorted by score in descending order
    Value  topScore = rootMoves[0].score;
    int    delta    = std::min(topScore - rootMoves[multiPV - 1].score, int(PawnValue));
    int    maxScore = -VALUE_INFINITE;
    double weakness = 120 - 2 * level;

    // Choose best move. For each move score we add two terms, both dependent on
    // weakness. One is deterministic and bigger for weaker levels, and one is
    // random. Then we choose the move with the resulting highest score.
    for (size_t i = 0; i < multiPV; ++i)
    {
        // This is our magic formula
        int push = (weakness * int(topScore - rootMoves[i].score)
                    + delta * (rng.rand<unsigned>() % int(weakness)))
                 / 128;

        if (rootMoves[i].score + push >= maxScore)
        {
            maxScore = rootMoves[i].score + push;
            best     = rootMoves[i].pv[0];
        }
    }

    return best;
}


// Used to print debug info and, more importantly, to detect
// when we are out of available time and thus stop the search.
void SearchManager::check_time(Search::Worker& worker) {
    if (--callsCnt > 0)
        return;

    // When using nodes, ensure checking rate is not lower than 0.1% of nodes
    callsCnt = worker.limits.nodes ? std::min(512, int(worker.limits.nodes / 1024)) : 512;

    static TimePoint lastInfoTime = now();

    TimePoint elapsed = tm.elapsed([&worker]() { return worker.threads.nodes_searched(); });
    TimePoint tick    = worker.limits.startTime + elapsed;

    if (tick - lastInfoTime >= 1000)
    {
        lastInfoTime = tick;
        dbg_print();
    }

    // We should not stop pondering until told so by the GUI
    if (ponder)
        return;

    if (
      // Later we rely on the fact that we can at least use the mainthread previous
      // root-search score and PV in a multithreaded environment to prove mated-in scores.
      worker.completedDepth >= 1
      && ((worker.limits.use_time_management() && (elapsed > tm.maximum() || stopOnPonderhit))
          || (worker.limits.movetime && elapsed >= worker.limits.movetime)
          || (worker.limits.nodes && worker.threads.nodes_searched() >= worker.limits.nodes)))
        worker.threads.stop = worker.threads.abortedSearch = true;
}

// Used to correct and extend PVs for moves that have a TB (but not a mate) score.
// Keeps the search based PV for as long as it is verified to maintain the game
// outcome, truncates afterwards. Finally, extends to mate the PV, providing a
// possible continuation (but not a proven mating line).
void syzygy_extend_pv(const OptionsMap&         options,
                      const Search::LimitsType& limits,
                      Position&                 pos,
                      RootMove&                 rootMove,
                      Value&                    v) {

    auto t_start      = std::chrono::steady_clock::now();
    int  moveOverhead = int(options["Move Overhead"]);
    bool rule50       = bool(options["Syzygy50MoveRule"]);

    // Do not use more than moveOverhead / 2 time, if time management is active
    auto time_abort = [&t_start, &moveOverhead, &limits]() -> bool {
        auto t_end = std::chrono::steady_clock::now();
        return limits.use_time_management()
            && 2 * std::chrono::duration<double, std::milli>(t_end - t_start).count()
                 > moveOverhead;
    };

    std::list<StateInfo> sts;

    // Step 0, do the rootMove, no correction allowed, as needed for MultiPV in TB.
    auto& stRoot = sts.emplace_back();
    pos.do_move(rootMove.pv[0], stRoot);
    int ply = 1;

    // Step 1, walk the PV to the last position in TB with correct decisive score
    while (size_t(ply) < rootMove.pv.size())
    {
        Move& pvMove = rootMove.pv[ply];

        RootMoves legalMoves;
        for (const auto& m : MoveList<LEGAL>(pos))
            legalMoves.emplace_back(m);

        Tablebases::Config config = Tablebases::rank_root_moves(options, pos, legalMoves);
        RootMove&          rm     = *std::find(legalMoves.begin(), legalMoves.end(), pvMove);

        if (legalMoves[0].tbRank != rm.tbRank)
            break;

        ply++;

        auto& st = sts.emplace_back();
        pos.do_move(pvMove, st);

        // Do not allow for repetitions or drawing moves along the PV in TB regime
        if (config.rootInTB && ((rule50 && pos.is_draw(ply)) || pos.is_repetition(ply)))
        {
            pos.undo_move(pvMove);
            ply--;
            break;
        }

        // Full PV shown will thus be validated and end in TB.
        // If we cannot validate the full PV in time, we do not show it.
        if (config.rootInTB && time_abort())
            break;
    }

    // Resize the PV to the correct part
    rootMove.pv.resize(ply);

    // Step 2, now extend the PV to mate, as if the user explored syzygy-tables.info
    // using top ranked moves (minimal DTZ), which gives optimal mates only for simple
    // endgames e.g. KRvK.
    while (!(rule50 && pos.is_draw(0)))
    {
        if (time_abort())
            break;

        RootMoves legalMoves;
        for (const auto& m : MoveList<LEGAL>(pos))
        {
            auto&     rm = legalMoves.emplace_back(m);
            StateInfo tmpSI;
            pos.do_move(m, tmpSI);
            // Give a score of each move to break DTZ ties restricting opponent mobility,
            // but not giving the opponent a capture.
            for (const auto& mOpp : MoveList<LEGAL>(pos))
                rm.tbRank -= pos.capture(mOpp) ? 100 : 1;
            pos.undo_move(m);
        }

        // Mate found
        if (legalMoves.size() == 0)
            break;

        // Sort moves according to their above assigned rank.
        // This will break ties for moves with equal DTZ in rank_root_moves.
        std::stable_sort(
          legalMoves.begin(), legalMoves.end(),
          [](const Search::RootMove& a, const Search::RootMove& b) { return a.tbRank > b.tbRank; });

        // The winning side tries to minimize DTZ, the losing side maximizes it
        Tablebases::Config config = Tablebases::rank_root_moves(options, pos, legalMoves, true);

        // If DTZ is not available we might not find a mate, so we bail out
        if (!config.rootInTB || config.cardinality > 0)
            break;

        ply++;

        Move& pvMove = legalMoves[0].pv[0];
        rootMove.pv.push_back(pvMove);
        auto& st = sts.emplace_back();
        pos.do_move(pvMove, st);
    }

    // Finding a draw in this function is an exceptional case, that cannot happen when rule50 is false or
    // during engine game play, since we have a winning score, and play correctly
    // with TB support. However, it can be that a position is draw due to the 50 move
    // rule if it has been been reached on the board with a non-optimal 50 move counter
    // (e.g. 8/8/6k1/3B4/3K4/4N3/8/8 w - - 54 106 ) which TB with dtz counter rounding
    // cannot always correctly rank. See also
    // https://github.com/official-stockfish/Stockfish/issues/5175#issuecomment-2058893495
    // We adjust the score to match the found PV. Note that a TB loss score can be
    // displayed if the engine did not find a drawing move yet, but eventually search
    // will figure it out (e.g. 1kq5/q2r4/5K2/8/8/8/8/7Q w - - 96 1 )
    if (pos.is_draw(0))
        v = VALUE_DRAW;

    // Undo the PV moves
    for (auto it = rootMove.pv.rbegin(); it != rootMove.pv.rend(); ++it)
        pos.undo_move(*it);

    // Inform if we couldn't get a full extension in time
    if (time_abort())
        sync_cout
          << "info string Syzygy based PV extension requires more time, increase Move Overhead as needed."
          << sync_endl;
}

void SearchManager::pv(Search::Worker&           worker,
                       const ThreadPool&         threads,
                       const TranspositionTable& tt,
                       Depth                     depth) {

    const auto nodes     = threads.nodes_searched();
    auto&      rootMoves = worker.rootMoves;
    auto&      pos       = worker.rootPos;
    size_t     pvIdx     = worker.pvIdx;
    size_t     multiPV   = std::min(size_t(worker.options["MultiPV"]), rootMoves.size());
    uint64_t   tbHits    = threads.tb_hits() + (worker.tbConfig.rootInTB ? rootMoves.size() : 0);

    for (size_t i = 0; i < multiPV; ++i)
    {
        bool updated = rootMoves[i].score != -VALUE_INFINITE;

        if (depth == 1 && !updated && i > 0)
            continue;

        Depth d = updated ? depth : std::max(1, depth - 1);
        Value v = updated ? rootMoves[i].uciScore : rootMoves[i].previousScore;

        if (v == -VALUE_INFINITE)
            v = VALUE_ZERO;

        bool tb = worker.tbConfig.rootInTB && std::abs(v) <= VALUE_TB;
        v       = tb ? rootMoves[i].tbScore : v;

        bool isExact = i != pvIdx || tb || !updated;  // tablebase- and previous-scores are exact

        // Potentially correct and extend the PV, and in exceptional cases v
        if (is_decisive(v) && std::abs(v) < VALUE_MATE_IN_MAX_PLY
            && ((!rootMoves[i].scoreLowerbound && !rootMoves[i].scoreUpperbound) || isExact))
            syzygy_extend_pv(worker.options, worker.limits, pos, rootMoves[i], v);

        std::string pv;
        for (Move m : rootMoves[i].pv)
            pv += UCIEngine::move(m, pos.is_chess960()) + " ";

        // Remove last whitespace
        if (!pv.empty())
            pv.pop_back();

        auto wdl   = worker.options["UCI_ShowWDL"] ? UCIEngine::wdl(v, pos) : "";
        auto bound = rootMoves[i].scoreLowerbound
                     ? "lowerbound"
                     : (rootMoves[i].scoreUpperbound ? "upperbound" : "");

        InfoFull info;

        info.depth    = d;
        info.selDepth = rootMoves[i].selDepth;
        info.multiPV  = i + 1;
        info.score    = {v, pos};
        info.wdl      = wdl;

        if (!isExact)
            info.bound = bound;

        TimePoint time = std::max(TimePoint(1), tm.elapsed_time());
        info.timeMs    = time;
        info.nodes     = nodes;
        info.nps       = nodes * 1000 / time;
        info.tbHits    = tbHits;
        info.pv        = pv;
        info.hashfull  = tt.hashfull();

        updates.onUpdateFull(info);
    }
}

// Called in case we have no ponder move before exiting the search,
// for instance, in case we stop the search during a fail high at root.
// We try hard to have a ponder move to return to the GUI,
// otherwise in case of 'ponder on' we have nothing to think about.
bool RootMove::extract_ponder_from_tt(const TranspositionTable& tt, Position& pos) {

    StateInfo st;
    ASSERT_ALIGNED(&st, Eval::NNUE::CacheLineSize);

    assert(pv.size() == 1);
    if (pv[0] == Move::none())
        return false;

    pos.do_move(pv[0], st, &tt);

    auto [ttHit, ttData, ttWriter] = tt.probe(pos.key());
    if (ttHit)
    {
        if (MoveList<LEGAL>(pos).contains(ttData.move))
            pv.push_back(ttData.move);
    }

    pos.undo_move(pv[0]);
    return pv.size() > 1;
}


}  // namespace Stockfish
