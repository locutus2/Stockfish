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

    return 7037 * pcv + 6671 * micv + 7631 * (wnpcv + bnpcv) + 6362 * cntcv;
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

    static constexpr int nonPawnWeight = 159;

    workerThread.pawnCorrectionHistory[pawn_structure_index<Correction>(pos)][us]
      << bonus * 104 / 128;
    workerThread.minorPieceCorrectionHistory[minor_piece_index(pos)][us] << bonus * 145 / 128;
    workerThread.nonPawnCorrectionHistory[WHITE][non_pawn_index<WHITE>(pos)][us]
      << bonus * nonPawnWeight / 128;
    workerThread.nonPawnCorrectionHistory[BLACK][non_pawn_index<BLACK>(pos)][us]
      << bonus * nonPawnWeight / 128;

    if (m.is_ok())
        (*(ss - 2)->continuationCorrectionHistory)[pos.piece_on(m.to_sq())][m.to_sq()]
          << bonus * 146 / 128;
}

// History and stats update bonus, based on depth
int stat_bonus(Depth d) { return std::min(154 * d - 102, 1661); }

// History and stats update malus, based on depth
int stat_malus(Depth d) { return std::min(831 * d - 269, 2666); }

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

    lowPlyHistory.fill(97);

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
            delta     = 5 + std::abs(rootMoves[pvIdx].meanSquaredScore) / 12991;
            Value avg = rootMoves[pvIdx].averageScore;
            alpha     = std::max(avg - delta, -VALUE_INFINITE);
            beta      = std::min(avg + delta, VALUE_INFINITE);

            // Adjust optimism based on root move's averageScore
            optimism[us]  = 141 * avg / (std::abs(avg) + 83);
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
    mainHistory.fill(63);
    lowPlyHistory.fill(108);
    captureHistory.fill(-631);
    pawnHistory.fill(-1210);
    pawnCorrectionHistory.fill(0);
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
                    h.fill(-479);

    for (size_t i = 1; i < reductions.size(); ++i)
        reductions[i] = int(2143 / 100.0 * std::log(i));

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
        && (cutNode == (ttData.value >= beta) || depth > 9))
    {
        // If ttMove is quiet, update move sorting heuristics on TT hit
        if (ttData.move && ttData.value >= beta)
        {
            // Bonus for a quiet ttMove that fails high
            if (!ttCapture)
                update_quiet_histories(pos, ss, *this, ttData.move, stat_bonus(depth) * 746 / 1024);

            // Extra penalty for early quiet moves of the previous ply
            if (prevSq != SQ_NONE && (ss - 1)->moveCount <= 2 && !priorCapture)
                update_continuation_histories(ss - 1, pos.piece_on(prevSq), prevSq,
                                              -stat_malus(depth + 1) * 1042 / 1024);
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
        int bonus = std::clamp(-10 * int((ss - 1)->staticEval + ss->staticEval), -1881, 1413) + 616;
        thisThread->mainHistory[~us][((ss - 1)->currentMove).from_to()] << bonus * 1151 / 1024;
        if (type_of(pos.piece_on(prevSq)) != PAWN && ((ss - 1)->currentMove).type_of() != PROMOTION)
            thisThread->pawnHistory[pawn_structure_index(pos)][pos.piece_on(prevSq)][prevSq]
              << bonus * 1107 / 1024;
    }

    // Set up the improving flag, which is true if current static evaluation is
    // bigger than the previous static evaluation at our turn (if we were in
    // check at our previous move we go back until we weren't in check) and is
    // false otherwise. The improving flag is used in various pruning heuristics.
    improving = ss->staticEval > (ss - 2)->staticEval;

    opponentWorsening = ss->staticEval + (ss - 1)->staticEval > 2;

    if (priorReduction >= 3 && !opponentWorsening)
        depth++;

    // Step 7. Razoring
    // If eval is really low, skip search entirely and return the qsearch value.
    // For PvNodes, we must have a guard against mates being returned.
    if (!PvNode && eval < alpha - 462 - 297 * depth * depth)
        return qsearch<NonPV>(pos, ss, alpha, beta);

    // Step 8. Futility pruning: child node
    // The depth condition is important for mate finding.
    if (!ss->ttPv && depth < 14
        && eval - futility_margin(depth, cutNode && !ss->ttHit, improving, opponentWorsening)
               - (ss - 1)->statScore / 310 + 40 - std::abs(correctionValue) / 131072
             >= beta
        && eval >= beta && (!ttData.move || ttCapture) && !is_loss(beta) && !is_win(eval))
        return beta + (eval - beta) / 3;

    // Step 9. Null move search with verification search
    if (cutNode && (ss - 1)->currentMove != Move::null() && eval >= beta
        && ss->staticEval >= beta - 20 * depth + 470 - 60 * improving && !excludedMove
        && pos.non_pawn_material(us) && ss->ply >= thisThread->nmpMinPly && !is_loss(beta))
    {
        assert(eval - beta >= 0);

        // Null move dynamic reduction based on depth and eval
        Depth R = std::min(int(eval - beta) / 215, 7) + depth / 3 + 5;

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
    if (((PvNode || cutNode) && depth >= 7 - 4 * PvNode) && !ttData.move)
        depth -= 2;

    // Step 11. ProbCut
    // If we have a good enough capture (or queen promotion) and a reduced search
    // returns a value much above beta, we can (almost) safely prune the previous move.
    probCutBeta = beta + 174 - 56 * improving;
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
    probCutBeta = beta + 412;
    if ((ttData.bound & BOUND_LOWER) && ttData.depth >= depth - 4 && ttData.value >= probCutBeta
        && !is_decisive(beta) && is_valid(ttData.value) && !is_decisive(ttData.value))
        return probCutBeta;

    const PieceToHistory* contHist[] = {
      (ss - 1)->continuationHistory, (ss - 2)->continuationHistory, (ss - 3)->continuationHistory,
      (ss - 4)->continuationHistory, (ss - 5)->continuationHistory, (ss - 6)->continuationHistory};


    MovePicker mp(pos, ttData.move, depth, &thisThread->mainHistory, &thisThread->lowPlyHistory,
                  &thisThread->captureHistory, contHist, &thisThread->pawnHistory, ss->ply);

    value = bestValue;

    int moveCount         = 0;
    int failedLMRResearch = 0;

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

        // Increase reduction for ttPv nodes (*Scaler)
        // Smaller or even negative value is better for short time controls
        // Bigger value is better for long time controls
        if (ss->ttPv)
            r += 1024;

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
                    Value futilityValue = ss->staticEval + 271 + 243 * lmrDepth
                                        + PieceValue[capturedPiece] + captHist / 7;
                    if (futilityValue <= alpha)
                        continue;
                }

                // SEE based pruning for captures and checks
                int seeHist = std::clamp(captHist / 37, -152 * depth, 141 * depth);
                if (!pos.see_ge(move, -156 * depth - seeHist))
                    continue;
            }
            else
            {
                int history =
                  (*contHist[0])[movedPiece][move.to_sq()]
                  + (*contHist[1])[movedPiece][move.to_sq()]
                  + thisThread->pawnHistory[pawn_structure_index(pos)][movedPiece][move.to_sq()];

                // Continuation history based pruning
                if (history < -3901 * depth)
                    continue;

                history += 2 * thisThread->mainHistory[us][move.from_to()];

                lmrDepth += history / 3459;

                Value futilityValue = ss->staticEval + (bestMove ? 47 : 137) + 142 * lmrDepth;

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
                if (!pos.see_ge(move, -25 * lmrDepth * lmrDepth))
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
                && depth >= 5 - (thisThread->completedDepth > 33) + ss->ttPv
                && is_valid(ttData.value) && !is_decisive(ttData.value)
                && (ttData.bound & BOUND_LOWER) && ttData.depth >= depth - 3)
            {
                Value singularBeta  = ttData.value - (52 + 74 * (ss->ttPv && !PvNode)) * depth / 64;
                Depth singularDepth = newDepth / 2;

                ss->excludedMove = move;
                value =
                  search<NonPV>(pos, ss, singularBeta - 1, singularBeta, singularDepth, cutNode);
                ss->excludedMove = Move::none();

                if (value < singularBeta)
                {
                    int corrValAdj   = std::abs(correctionValue) / 262144;
                    int doubleMargin = 249 * PvNode - 194 * !ttCapture - corrValAdj;
                    int tripleMargin =
                      94 + 287 * PvNode - 249 * !ttCapture + 99 * ss->ttPv - corrValAdj;

                    extension = 1 + (value < singularBeta - doubleMargin)
                              + (value < singularBeta - tripleMargin);

                    depth += (depth < 15);
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

            // Extension for capturing the previous moved piece
            else if (PvNode && move.to_sq() == prevSq
                     && thisThread->captureHistory[movedPiece][move.to_sq()]
                                                  [type_of(pos.piece_on(move.to_sq()))]
                          > 4126)
                extension = 1;
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
            r -= 2061 + (ttData.value > alpha) * 965 + (ttData.depth >= depth) * 960;

        if (PvNode)
            r -= 1018;

        // These reduction adjustments have no proven non-linear scaling

        r += 307 - moveCount * 64;

        r -= std::abs(correctionValue) / 34112;

        // Increase reduction for cut nodes
        if (cutNode)
            r += 2355 - (ttData.depth >= depth && ss->ttPv) * 1141;

        // Increase reduction if ttMove is a capture but the current move is not a capture
        if (ttCapture && !capture)
            r += 1087 + (depth < 8) * 990;

        // Increase reduction if next ply has a lot of fail high
        if ((ss + 1)->cutoffCnt > 3)
            r += 940 + allNode * 887;

        // For first picked move (ttMove) reduce reduction
        else if (move == ttData.move)
            r -= 1960;

        if (capture)
            ss->statScore =
              7 * int(PieceValue[pos.captured_piece()])
              + thisThread->captureHistory[movedPiece][move.to_sq()][type_of(pos.captured_piece())]
              - 4666;
        else
            ss->statScore = 2 * thisThread->mainHistory[us][move.from_to()]
                          + (*contHist[0])[movedPiece][move.to_sq()]
                          + (*contHist[1])[movedPiece][move.to_sq()] - 3874;

        // Decrease/increase reduction for moves with a good/bad history
        r -= ss->statScore * 1451 / 16384;

        //r += failedLMRResearch * 1024;

        // Step 17. Late moves reduction / extension (LMR)
        if (depth >= 2 && moveCount > 1)
        {
            // In general we want to cap the LMR depth search at newDepth, but when
            // reduction is negative, we allow this move a limited search extension
            // beyond the first move depth.
            // To prevent problems when the max value is less than the min value,
            // std::clamp has been replaced by a more robust implementation.

            bool C = ss->ttHit;

            Depth d = std::max(
              1, std::min(newDepth - r / 1024, newDepth + !allNode + (PvNode && !bestMove)));

            ss->reduction = newDepth - d;

            value         = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, d, true);
            ss->reduction = 0;


            // Do a full-depth search when reduced LMR search fails high
            if (value > alpha && d < newDepth)
            {
                // Adjust full-depth search based on LMR results - if the result was
                // good enough search deeper, if it was bad enough search shallower.
                const bool doDeeperSearch    = value > (bestValue + 40 + 2 * newDepth);
                const bool doShallowerSearch = value < bestValue + 10;

                newDepth += doDeeperSearch - doShallowerSearch;

                if (newDepth > d)
                {
                    value = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode);

                    dbg_hit_on(value <= alpha, failedLMRResearch + 1000*C);
                    if (value <= alpha)
                        failedLMRResearch++;

                    /*
                    bool C = ss->ttHit;
                    Hit #0: Total 961734 Hits 197528 Hit Rate (%) 20.5387
                    Hit #1: Total 68623 Hits 26538 Hit Rate (%) 38.6722
                    Hit #2: Total 10981 Hits 5896 Hit Rate (%) 53.6927
                    Hit #3: Total 2505 Hits 1599 Hit Rate (%) 63.8323
                    Hit #4: Total 736 Hits 508 Hit Rate (%) 69.0217
                    Hit #5: Total 216 Hits 143 Hit Rate (%) 66.2037
                    Hit #6: Total 59 Hits 48 Hit Rate (%) 81.3559
                    Hit #7: Total 15 Hits 13 Hit Rate (%) 86.6667
                    Hit #8: Total 2 Hits 2 Hit Rate (%) 100
                    Hit #9: Total 2 Hits 2 Hit Rate (%) 100
                    Hit #1000: Total 3577398 Hits 888196 Hit Rate (%) 24.828
                    Hit #1001: Total 276320 Hits 119958 Hit Rate (%) 43.4127
                    Hit #1002: Total 43299 Hits 24825 Hit Rate (%) 57.3339
                    Hit #1003: Total 9627 Hits 6329 Hit Rate (%) 65.7422
                    Hit #1004: Total 2528 Hits 1803 Hit Rate (%) 71.3212
                    Hit #1005: Total 759 Hits 580 Hit Rate (%) 76.4163
                    Hit #1006: Total 237 Hits 185 Hit Rate (%) 78.0591
                    Hit #1007: Total 76 Hits 63 Hit Rate (%) 82.8947
                    Hit #1008: Total 30 Hits 27 Hit Rate (%) 90
                    Hit #1009: Total 9 Hits 6 Hit Rate (%) 66.6667
                    Hit #1010: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #1011: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1012: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1013: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = ss->statScore > 0;
                    Hit #0: Total 2591040 Hits 750678 Hit Rate (%) 28.9721
                    Hit #1: Total 265551 Hits 121827 Hit Rate (%) 45.8771
                    Hit #2: Total 47341 Hits 27734 Hit Rate (%) 58.5835
                    Hit #3: Total 11111 Hits 7372 Hit Rate (%) 66.3487
                    Hit #4: Total 3093 Hits 2218 Hit Rate (%) 71.7103
                    Hit #5: Total 936 Hits 699 Hit Rate (%) 74.6795
                    Hit #6: Total 291 Hits 231 Hit Rate (%) 79.3814
                    Hit #7: Total 87 Hits 74 Hit Rate (%) 85.0575
                    Hit #8: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 1948092 Hits 335046 Hit Rate (%) 17.1987
                    Hit #1001: Total 79392 Hits 24669 Hit Rate (%) 31.0724
                    Hit #1002: Total 6939 Hits 2987 Hit Rate (%) 43.0465
                    Hit #1003: Total 1021 Hits 556 Hit Rate (%) 54.4564
                    Hit #1004: Total 171 Hits 93 Hit Rate (%) 54.386
                    Hit #1005: Total 39 Hits 24 Hit Rate (%) 61.5385
                    Hit #1006: Total 5 Hits 2 Hit Rate (%) 40
                    Hit #1007: Total 4 Hits 2 Hit Rate (%) 50

                    bool C = capture && givesCheck;
                    Hit #0: Total 4489044 Hits 1079200 Hit Rate (%) 24.0408
                    Hit #1: Total 343242 Hits 146081 Hit Rate (%) 42.5592
                    Hit #2: Total 54131 Hits 30674 Hit Rate (%) 56.6662
                    Hit #3: Total 12109 Hits 7918 Hit Rate (%) 65.3894
                    Hit #4: Total 3261 Hits 2309 Hit Rate (%) 70.8065
                    Hit #5: Total 974 Hits 722 Hit Rate (%) 74.1273
                    Hit #6: Total 296 Hits 233 Hit Rate (%) 78.7162
                    Hit #7: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #8: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 50088 Hits 6524 Hit Rate (%) 13.0251
                    Hit #1001: Total 1701 Hits 415 Hit Rate (%) 24.3974
                    Hit #1002: Total 149 Hits 47 Hit Rate (%) 31.5436
                    Hit #1003: Total 23 Hits 10 Hit Rate (%) 43.4783
                    Hit #1004: Total 3 Hits 2 Hit Rate (%) 66.6667
                    Hit #1005: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = PvNode && capture;
                    Hit #0: Total 4538989 Hits 1085705 Hit Rate (%) 23.9195
                    Hit #1: Total 344943 Hits 146496 Hit Rate (%) 42.4696
                    Hit #2: Total 54280 Hits 30721 Hit Rate (%) 56.5973
                    Hit #3: Total 12132 Hits 7928 Hit Rate (%) 65.3478
                    Hit #4: Total 3264 Hits 2311 Hit Rate (%) 70.8027
                    Hit #5: Total 975 Hits 723 Hit Rate (%) 74.1538
                    Hit #6: Total 296 Hits 233 Hit Rate (%) 78.7162
                    Hit #7: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #8: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 143 Hits 19 Hit Rate (%) 13.2867

                    bool C = (ss+1)->cutoffCnt > 0;
                    Hit #0: Total 1072536 Hits 205696 Hit Rate (%) 19.1785
                    Hit #1: Total 143 Hits 79 Hit Rate (%) 55.2448
                    Hit #2: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #1000: Total 3466596 Hits 880028 Hit Rate (%) 25.3859
                    Hit #1001: Total 344800 Hits 146417 Hit Rate (%) 42.4643
                    Hit #1002: Total 54276 Hits 30718 Hit Rate (%) 56.5959
                    Hit #1003: Total 12132 Hits 7928 Hit Rate (%) 65.3478
                    Hit #1004: Total 3264 Hits 2311 Hit Rate (%) 70.8027
                    Hit #1005: Total 975 Hits 723 Hit Rate (%) 74.1538
                    Hit #1006: Total 296 Hits 233 Hit Rate (%) 78.7162
                    Hit #1007: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #1008: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #1009: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #1010: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #1011: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1012: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1013: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = (ss-1)->currentMove == Move::null();
                    Hit #0: Total 4521750 Hits 1081632 Hit Rate (%) 23.9207
                    Hit #1: Total 344319 Hits 146228 Hit Rate (%) 42.4688
                    Hit #2: Total 54239 Hits 30697 Hit Rate (%) 56.5958
                    Hit #3: Total 12123 Hits 7921 Hit Rate (%) 65.3386
                    Hit #4: Total 3264 Hits 2311 Hit Rate (%) 70.8027
                    Hit #5: Total 975 Hits 723 Hit Rate (%) 74.1538
                    Hit #6: Total 296 Hits 233 Hit Rate (%) 78.7162
                    Hit #7: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #8: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 17382 Hits 4092 Hit Rate (%) 23.5416
                    Hit #1001: Total 624 Hits 268 Hit Rate (%) 42.9487
                    Hit #1002: Total 41 Hits 24 Hit Rate (%) 58.5366
                    Hit #1003: Total 9 Hits 7 Hit Rate (%) 77.7778

                    bool C = (ss+1)->cutoffCnt > 1;
                    Hit #0: Total 2676119 Hits 565757 Hit Rate (%) 21.141
                    Hit #1: Total 42480 Hits 14449 Hit Rate (%) 34.0137
                    Hit #2: Total 37 Hits 25 Hit Rate (%) 67.5676
                    Hit #3: Total 3 Hits 3 Hit Rate (%) 100
                    Hit #1000: Total 1863013 Hits 519967 Hit Rate (%) 27.91
                    Hit #1001: Total 302463 Hits 132047 Hit Rate (%) 43.6572
                    Hit #1002: Total 54243 Hits 30696 Hit Rate (%) 56.5898
                    Hit #1003: Total 12129 Hits 7925 Hit Rate (%) 65.3393
                    Hit #1004: Total 3264 Hits 2311 Hit Rate (%) 70.8027
                    Hit #1005: Total 975 Hits 723 Hit Rate (%) 74.1538
                    Hit #1006: Total 296 Hits 233 Hit Rate (%) 78.7162
                    Hit #1007: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #1008: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #1009: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #1010: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #1011: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1012: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1013: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = capture || givesCheck;
                    Hit #0: Total 3520003 Hits 876136 Hit Rate (%) 24.8902
                    Hit #1: Total 297514 Hits 129043 Hit Rate (%) 43.3738
                    Hit #2: Total 49017 Hits 28245 Hit Rate (%) 57.6229
                    Hit #3: Total 11261 Hits 7470 Hit Rate (%) 66.3351
                    Hit #4: Total 3103 Hits 2217 Hit Rate (%) 71.447
                    Hit #5: Total 921 Hits 692 Hit Rate (%) 75.1357
                    Hit #6: Total 284 Hits 225 Hit Rate (%) 79.2254
                    Hit #7: Total 89 Hits 74 Hit Rate (%) 83.1461
                    Hit #8: Total 31 Hits 28 Hit Rate (%) 90.3226
                    Hit #9: Total 10 Hits 7 Hit Rate (%) 70
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 1019129 Hits 209588 Hit Rate (%) 20.5654
                    Hit #1001: Total 47429 Hits 17453 Hit Rate (%) 36.7982
                    Hit #1002: Total 5263 Hits 2476 Hit Rate (%) 47.0454
                    Hit #1003: Total 871 Hits 458 Hit Rate (%) 52.5832
                    Hit #1004: Total 161 Hits 94 Hit Rate (%) 58.3851
                    Hit #1005: Total 54 Hits 31 Hit Rate (%) 57.4074
                    Hit #1006: Total 12 Hits 8 Hit Rate (%) 66.6667
                    Hit #1007: Total 2 Hits 2 Hit Rate (%) 100
                    Hit #1008: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1009: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = (ss+1)->cutoffCnt > 2;
                    Hit #0: Total 3378637 Hits 745013 Hit Rate (%) 22.0507
                    Hit #1: Total 137342 Hits 51288 Hit Rate (%) 37.3433
                    Hit #2: Total 3866 Hits 1824 Hit Rate (%) 47.1805
                    Hit #3: Total 15 Hits 11 Hit Rate (%) 73.3333
                    Hit #4: Total 2 Hits 2 Hit Rate (%) 100
                    Hit #1000: Total 1160495 Hits 340711 Hit Rate (%) 29.3591
                    Hit #1001: Total 207601 Hits 95208 Hit Rate (%) 45.8611
                    Hit #1002: Total 50414 Hits 28897 Hit Rate (%) 57.3194
                    Hit #1003: Total 12117 Hits 7917 Hit Rate (%) 65.338
                    Hit #1004: Total 3262 Hits 2309 Hit Rate (%) 70.7848
                    Hit #1005: Total 975 Hits 723 Hit Rate (%) 74.1538
                    Hit #1006: Total 296 Hits 233 Hit Rate (%) 78.7162
                    Hit #1007: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #1008: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #1009: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #1010: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #1011: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1012: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1013: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = (ss+1)->cutoffCnt > 3;
                    Hit #0: Total 3748208 Hits 845539 Hit Rate (%) 22.5585
                    Hit #1: Total 200565 Hits 77806 Hit Rate (%) 38.7934
                    Hit #2: Total 15434 Hits 7978 Hit Rate (%) 51.6911
                    Hit #3: Total 559 Hits 338 Hit Rate (%) 60.4651
                    Hit #4: Total 5 Hits 5 Hit Rate (%) 100
                    Hit #1000: Total 790924 Hits 240185 Hit Rate (%) 30.3676
                    Hit #1001: Total 144378 Hits 68690 Hit Rate (%) 47.5765
                    Hit #1002: Total 38846 Hits 22743 Hit Rate (%) 58.5466
                    Hit #1003: Total 11573 Hits 7590 Hit Rate (%) 65.5837
                    Hit #1004: Total 3259 Hits 2306 Hit Rate (%) 70.7579
                    Hit #1005: Total 975 Hits 723 Hit Rate (%) 74.1538
                    Hit #1006: Total 296 Hits 233 Hit Rate (%) 78.7162
                    Hit #1007: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #1008: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #1009: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #1010: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #1011: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1012: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1013: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = PvNode || capture;
                    Hit #0: Total 4084601 Hits 1011253 Hit Rate (%) 24.7577
                    Hit #1: Total 330492 Hits 142272 Hit Rate (%) 43.0485
                    Hit #2: Total 52870 Hits 30223 Hit Rate (%) 57.1647
                    Hit #3: Total 11934 Hits 7849 Hit Rate (%) 65.7701
                    Hit #4: Total 3231 Hits 2295 Hit Rate (%) 71.0306
                    Hit #5: Total 972 Hits 722 Hit Rate (%) 74.2798
                    Hit #6: Total 295 Hits 233 Hit Rate (%) 78.9831
                    Hit #7: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #8: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 454531 Hits 74471 Hit Rate (%) 16.3841
                    Hit #1001: Total 14451 Hits 4224 Hit Rate (%) 29.2298
                    Hit #1002: Total 1410 Hits 498 Hit Rate (%) 35.3191
                    Hit #1003: Total 198 Hits 79 Hit Rate (%) 39.899
                    Hit #1004: Total 33 Hits 16 Hit Rate (%) 48.4848
                    Hit #1005: Total 3 Hits 1 Hit Rate (%) 33.3333
                    Hit #1006: Total 1 Hits 0 Hit Rate (%) 0

                    bool C = ss->ttPv;
                    Hit #0: Total 4130588 Hits 952062 Hit Rate (%) 23.0491
                    Hit #1: Total 315043 Hits 131435 Hit Rate (%) 41.7197
                    Hit #2: Total 50305 Hits 28223 Hit Rate (%) 56.1038
                    Hit #3: Total 11436 Hits 7458 Hit Rate (%) 65.2151
                    Hit #4: Total 3122 Hits 2190 Hit Rate (%) 70.1473
                    Hit #5: Total 932 Hits 687 Hit Rate (%) 73.7124
                    Hit #6: Total 283 Hits 223 Hit Rate (%) 78.7986
                    Hit #7: Total 89 Hits 74 Hit Rate (%) 83.1461
                    Hit #8: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 408544 Hits 133662 Hit Rate (%) 32.7167
                    Hit #1001: Total 29900 Hits 15061 Hit Rate (%) 50.3712
                    Hit #1002: Total 3975 Hits 2498 Hit Rate (%) 62.8428
                    Hit #1003: Total 696 Hits 470 Hit Rate (%) 67.5287
                    Hit #1004: Total 142 Hits 121 Hit Rate (%) 85.2113
                    Hit #1005: Total 43 Hits 36 Hit Rate (%) 83.7209
                    Hit #1006: Total 13 Hits 10 Hit Rate (%) 76.9231
                    Hit #1007: Total 2 Hits 2 Hit Rate (%) 100

                    bool C = ttCapture;
                    Hit #0: Total 4227753 Hits 1025377 Hit Rate (%) 24.2535
                    Hit #1: Total 330030 Hits 140579 Hit Rate (%) 42.5958
                    Hit #2: Total 52537 Hits 29752 Hit Rate (%) 56.6306
                    Hit #3: Total 11758 Hits 7673 Hit Rate (%) 65.2577
                    Hit #4: Total 3161 Hits 2231 Hit Rate (%) 70.5789
                    Hit #5: Total 938 Hits 691 Hit Rate (%) 73.6674
                    Hit #6: Total 284 Hits 223 Hit Rate (%) 78.5211
                    Hit #7: Total 86 Hits 72 Hit Rate (%) 83.7209
                    Hit #8: Total 30 Hits 27 Hit Rate (%) 90
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 311379 Hits 60347 Hit Rate (%) 19.3806
                    Hit #1001: Total 14913 Hits 5917 Hit Rate (%) 39.6768
                    Hit #1002: Total 1743 Hits 969 Hit Rate (%) 55.5938
                    Hit #1003: Total 374 Hits 255 Hit Rate (%) 68.1818
                    Hit #1004: Total 103 Hits 80 Hit Rate (%) 77.6699
                    Hit #1005: Total 37 Hits 32 Hit Rate (%) 86.4865
                    Hit #1006: Total 12 Hits 10 Hit Rate (%) 83.3333
                    Hit #1007: Total 5 Hits 4 Hit Rate (%) 80
                    Hit #1008: Total 2 Hits 2 Hit Rate (%) 100

                    bool C = givesCheck;
                    Hit #0: Total 3906913 Hits 939797 Hit Rate (%) 24.0547
                    Hit #1: Total 309567 Hits 132612 Hit Rate (%) 42.8379
                    Hit #2: Total 50220 Hits 28673 Hit Rate (%) 57.0948
                    Hit #3: Total 11429 Hits 7534 Hit Rate (%) 65.92
                    Hit #4: Total 3132 Hits 2230 Hit Rate (%) 71.2005
                    Hit #5: Total 923 Hits 692 Hit Rate (%) 74.9729
                    Hit #6: Total 285 Hits 225 Hit Rate (%) 78.9474
                    Hit #7: Total 89 Hits 74 Hit Rate (%) 83.1461
                    Hit #8: Total 31 Hits 28 Hit Rate (%) 90.3226
                    Hit #9: Total 10 Hits 7 Hit Rate (%) 70
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 632219 Hits 145927 Hit Rate (%) 23.0817
                    Hit #1001: Total 35376 Hits 13884 Hit Rate (%) 39.2469
                    Hit #1002: Total 4060 Hits 2048 Hit Rate (%) 50.4433
                    Hit #1003: Total 703 Hits 394 Hit Rate (%) 56.0455
                    Hit #1004: Total 132 Hits 81 Hit Rate (%) 61.3636
                    Hit #1005: Total 52 Hits 31 Hit Rate (%) 59.6154
                    Hit #1006: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #1007: Total 2 Hits 2 Hit Rate (%) 100
                    Hit #1008: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1009: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = ss->inCheck;
                    Hit #0: Total 4288121 Hits 1035212 Hit Rate (%) 24.1414
                    Hit #1: Total 336618 Hits 143128 Hit Rate (%) 42.5194
                    Hit #2: Total 53675 Hits 30358 Hit Rate (%) 56.5589
                    Hit #3: Total 12047 Hits 7857 Hit Rate (%) 65.2196
                    Hit #4: Total 3247 Hits 2296 Hit Rate (%) 70.7114
                    Hit #5: Total 974 Hits 722 Hit Rate (%) 74.1273
                    Hit #6: Total 296 Hits 233 Hit Rate (%) 78.7162
                    Hit #7: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #8: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 251011 Hits 50512 Hit Rate (%) 20.1234
                    Hit #1001: Total 8325 Hits 3368 Hit Rate (%) 40.4565
                    Hit #1002: Total 605 Hits 363 Hit Rate (%) 60
                    Hit #1003: Total 85 Hits 71 Hit Rate (%) 83.5294
                    Hit #1004: Total 17 Hits 15 Hit Rate (%) 88.2353
                    Hit #1005: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = capture;
                    Hit #0: Total 4102134 Hits 1015539 Hit Rate (%) 24.7564
                    Hit #1: Total 331189 Hits 142512 Hit Rate (%) 43.0304
                    Hit #2: Total 52928 Hits 30246 Hit Rate (%) 57.1456
                    Hit #3: Total 11941 Hits 7854 Hit Rate (%) 65.7734
                    Hit #4: Total 3232 Hits 2296 Hit Rate (%) 71.0396
                    Hit #5: Total 972 Hits 722 Hit Rate (%) 74.2798
                    Hit #6: Total 295 Hits 233 Hit Rate (%) 78.9831
                    Hit #7: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #8: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 436998 Hits 70185 Hit Rate (%) 16.0607
                    Hit #1001: Total 13754 Hits 3984 Hit Rate (%) 28.9661
                    Hit #1002: Total 1352 Hits 475 Hit Rate (%) 35.1331
                    Hit #1003: Total 191 Hits 74 Hit Rate (%) 38.7435
                    Hit #1004: Total 32 Hits 15 Hit Rate (%) 46.875
                    Hit #1005: Total 3 Hits 1 Hit Rate (%) 33.3333
                    Hit #1006: Total 1 Hits 0 Hit Rate (%) 0
                     
                    bool C = priorCapture;
                    Hit #0: Total 3938495 Hits 952915 Hit Rate (%) 24.1949
                    Hit #1: Total 311586 Hits 132466 Hit Rate (%) 42.5135
                    Hit #2: Total 49282 Hits 27889 Hit Rate (%) 56.5906
                    Hit #3: Total 10989 Hits 7179 Hit Rate (%) 65.329
                    Hit #4: Total 2961 Hits 2106 Hit Rate (%) 71.1246
                    Hit #5: Total 892 Hits 659 Hit Rate (%) 73.8789
                    Hit #6: Total 274 Hits 219 Hit Rate (%) 79.927
                    Hit #7: Total 86 Hits 72 Hit Rate (%) 83.7209
                    Hit #8: Total 30 Hits 27 Hit Rate (%) 90
                    Hit #9: Total 10 Hits 7 Hit Rate (%) 70
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 600637 Hits 132809 Hit Rate (%) 22.1114
                    Hit #1001: Total 33357 Hits 14030 Hit Rate (%) 42.0601
                    Hit #1002: Total 4998 Hits 2832 Hit Rate (%) 56.6627
                    Hit #1003: Total 1143 Hits 749 Hit Rate (%) 65.5293
                    Hit #1004: Total 303 Hits 205 Hit Rate (%) 67.6568
                    Hit #1005: Total 83 Hits 64 Hit Rate (%) 77.1084
                    Hit #1006: Total 22 Hits 14 Hit Rate (%) 63.6364
                    Hit #1007: Total 5 Hits 4 Hit Rate (%) 80
                    Hit #1008: Total 2 Hits 2 Hit Rate (%) 100
                    Hit #1009: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = improving;
                    Hit #0: Total 1692006 Hits 380435 Hit Rate (%) 22.4843
                    Hit #1: Total 66624 Hits 26916 Hit Rate (%) 40.3999
                    Hit #2: Total 5897 Hits 3390 Hit Rate (%) 57.4869
                    Hit #3: Total 874 Hits 610 Hit Rate (%) 69.7941
                    Hit #4: Total 167 Hits 140 Hit Rate (%) 83.8323
                    Hit #5: Total 49 Hits 39 Hit Rate (%) 79.5918
                    Hit #6: Total 10 Hits 10 Hit Rate (%) 100
                    Hit #7: Total 2 Hits 2 Hit Rate (%) 100
                    Hit #1000: Total 2847126 Hits 705289 Hit Rate (%) 24.772
                    Hit #1001: Total 278319 Hits 119580 Hit Rate (%) 42.9651
                    Hit #1002: Total 48383 Hits 27331 Hit Rate (%) 56.4888
                    Hit #1003: Total 11258 Hits 7318 Hit Rate (%) 65.0027
                    Hit #1004: Total 3097 Hits 2171 Hit Rate (%) 70.1001
                    Hit #1005: Total 926 Hits 684 Hit Rate (%) 73.8661
                    Hit #1006: Total 286 Hits 223 Hit Rate (%) 77.972
                    Hit #1007: Total 89 Hits 74 Hit Rate (%) 83.1461
                    Hit #1008: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #1009: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #1010: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #1011: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1012: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1013: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = allNode;
                    Hit #0: Total 3276137 Hits 777689 Hit Rate (%) 23.738
                    Hit #1: Total 294760 Hits 126668 Hit Rate (%) 42.9733
                    Hit #2: Total 49938 Hits 28293 Hit Rate (%) 56.6563
                    Hit #3: Total 11444 Hits 7438 Hit Rate (%) 64.9948
                    Hit #4: Total 3090 Hits 2175 Hit Rate (%) 70.3883
                    Hit #5: Total 924 Hits 687 Hit Rate (%) 74.3506
                    Hit #6: Total 283 Hits 221 Hit Rate (%) 78.0919
                    Hit #7: Total 88 Hits 73 Hit Rate (%) 82.9545
                    Hit #8: Total 31 Hits 28 Hit Rate (%) 90.3226
                    Hit #9: Total 10 Hits 7 Hit Rate (%) 70
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 1262995 Hits 308035 Hit Rate (%) 24.3892
                    Hit #1001: Total 50183 Hits 19828 Hit Rate (%) 39.5114
                    Hit #1002: Total 4342 Hits 2428 Hit Rate (%) 55.9189
                    Hit #1003: Total 688 Hits 490 Hit Rate (%) 71.2209
                    Hit #1004: Total 174 Hits 136 Hit Rate (%) 78.1609
                    Hit #1005: Total 51 Hits 36 Hit Rate (%) 70.5882
                    Hit #1006: Total 13 Hits 12 Hit Rate (%) 92.3077
                    Hit #1007: Total 3 Hits 3 Hit Rate (%) 100
                    Hit #1008: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1009: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = cutNode;
                    Hit #0: Total 1280671 Hits 312340 Hit Rate (%) 24.3888
                    Hit #1: Total 50880 Hits 20068 Hit Rate (%) 39.4418
                    Hit #2: Total 4400 Hits 2451 Hit Rate (%) 55.7045
                    Hit #3: Total 695 Hits 495 Hit Rate (%) 71.223
                    Hit #4: Total 175 Hits 137 Hit Rate (%) 78.2857
                    Hit #5: Total 51 Hits 36 Hit Rate (%) 70.5882
                    Hit #6: Total 13 Hits 12 Hit Rate (%) 92.3077
                    Hit #7: Total 3 Hits 3 Hit Rate (%) 100
                    Hit #8: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #9: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 3258461 Hits 773384 Hit Rate (%) 23.7346
                    Hit #1001: Total 294063 Hits 126428 Hit Rate (%) 42.9935
                    Hit #1002: Total 49880 Hits 28270 Hit Rate (%) 56.676
                    Hit #1003: Total 11437 Hits 7433 Hit Rate (%) 64.9908
                    Hit #1004: Total 3089 Hits 2174 Hit Rate (%) 70.3788
                    Hit #1005: Total 924 Hits 687 Hit Rate (%) 74.3506
                    Hit #1006: Total 283 Hits 221 Hit Rate (%) 78.0919
                    Hit #1007: Total 88 Hits 73 Hit Rate (%) 82.9545
                    Hit #1008: Total 31 Hits 28 Hit Rate (%) 90.3226
                    Hit #1009: Total 10 Hits 7 Hit Rate (%) 70
                    Hit #1010: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #1011: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1012: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1013: Total 1 Hits 1 Hit Rate (%) 100

                    bool C = PvNode;
                    Hit #0: Total 4521456 Hits 1081419 Hit Rate (%) 23.9175
                    Hit #1: Total 344246 Hits 146256 Hit Rate (%) 42.4859
                    Hit #2: Total 54222 Hits 30698 Hit Rate (%) 56.6154
                    Hit #3: Total 12125 Hits 7923 Hit Rate (%) 65.3443
                    Hit #4: Total 3263 Hits 2310 Hit Rate (%) 70.7937
                    Hit #5: Total 975 Hits 723 Hit Rate (%) 74.1538
                    Hit #6: Total 296 Hits 233 Hit Rate (%) 78.7162
                    Hit #7: Total 91 Hits 76 Hit Rate (%) 83.5165
                    Hit #8: Total 32 Hits 29 Hit Rate (%) 90.625
                    Hit #9: Total 11 Hits 8 Hit Rate (%) 72.7273
                    Hit #10: Total 4 Hits 3 Hit Rate (%) 75
                    Hit #11: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #12: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #13: Total 1 Hits 1 Hit Rate (%) 100
                    Hit #1000: Total 17676 Hits 4305 Hit Rate (%) 24.3551
                    Hit #1001: Total 697 Hits 240 Hit Rate (%) 34.4333
                    Hit #1002: Total 58 Hits 23 Hit Rate (%) 39.6552
                    Hit #1003: Total 7 Hits 5 Hit Rate (%) 71.4286
                    Hit #1004: Total 1 Hits 1 Hit Rate (%) 100
                     * */
                }

                // Post LMR continuation history updates
                int bonus = (value >= beta) * 2048;
                update_continuation_histories(ss, movedPiece, move.to_sq(), bonus);
            }
        }

        // Step 18. Full-depth search when LMR is skipped
        else if (!PvNode || moveCount > 1)
        {
            // Increase reduction if ttMove is not present
            if (!ttData.move)
                r += 2111;

            // Note that if expected reduction is high, we reduce search depth here
            value =
              -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, newDepth - (r > 3444), !cutNode);
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
                    ss->cutoffCnt += (extension < 2);
                    assert(value >= beta);  // Fail high
                    break;
                }
                else
                {
                    // Reduce other moves if we have found at least one score improvement
                    if (depth > 2 && depth < 14 && !is_decisive(value))
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
        int bonusScale = (118 * (depth > 5) + 37 * !allNode + 169 * ((ss - 1)->moveCount > 8)
                          + 128 * (!ss->inCheck && bestValue <= ss->staticEval - 102)
                          + 115 * (!(ss - 1)->inCheck && bestValue <= -(ss - 1)->staticEval - 82)
                          + 80 * ((ss - 1)->isTTMove) + std::min(-(ss - 1)->statScore / 106, 318));

        bonusScale = std::max(bonusScale, 0);

        const int scaledBonus = stat_bonus(depth) * bonusScale;

        update_continuation_histories(ss - 1, pos.piece_on(prevSq), prevSq,
                                      scaledBonus * 436 / 32768);

        thisThread->mainHistory[~us][((ss - 1)->currentMove).from_to()]
          << scaledBonus * 207 / 32768;

        if (type_of(pos.piece_on(prevSq)) != PAWN && ((ss - 1)->currentMove).type_of() != PROMOTION)
            thisThread->pawnHistory[pawn_structure_index(pos)][pos.piece_on(prevSq)][prevSq]
              << scaledBonus * 1195 / 32768;
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
    Color us = pos.side_to_move();

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

        futilityBase = ss->staticEval + 301;
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
        if (!is_loss(bestValue) && pos.non_pawn_material(us))
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
                     <= 5228)
                continue;

            // Do not search moves with bad enough SEE values
            if (!pos.see_ge(move, -80))
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
    return reductionScale - delta * 768 / rootDelta + !i * reductionScale * 108 / 300 + 1168;
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

    int bonus = stat_bonus(depth) + 300 * isTTMove;
    int malus = stat_malus(depth) - 34 * (moveCount - 1);

    if (!pos.capture_stage(bestMove))
    {
        update_quiet_histories(pos, ss, workerThread, bestMove, bonus * 1216 / 1024);

        // Decrease stats for all non-best quiet moves
        for (Move move : quietsSearched)
            update_quiet_histories(pos, ss, workerThread, move, -malus * 1062 / 1024);
    }
    else
    {
        // Increase stats for the best move in case it was a capture move
        captured = type_of(pos.piece_on(bestMove.to_sq()));
        captureHistory[moved_piece][bestMove.to_sq()][captured] << bonus * 1272 / 1024;
    }

    // Extra penalty for a quiet early move that was not a TT move in
    // previous ply when it gets refuted.
    if (prevSq != SQ_NONE && ((ss - 1)->moveCount == 1 + (ss - 1)->ttHit) && !pos.captured_piece())
        update_continuation_histories(ss - 1, pos.piece_on(prevSq), prevSq, -malus * 966 / 1024);

    // Decrease stats for all non-best capture moves
    for (Move move : capturesSearched)
    {
        moved_piece = pos.moved_piece(move);
        captured    = type_of(pos.piece_on(move.to_sq()));
        captureHistory[moved_piece][move.to_sq()][captured] << -malus * 1205 / 1024;
    }
}


// Updates histories of the move pairs formed by moves
// at ply -1, -2, -3, -4, and -6 with current move.
void update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus) {
    static constexpr std::array<ConthistBonus, 6> conthist_bonuses = {
      {{1, 1025}, {2, 621}, {3, 325}, {4, 512}, {5, 122}, {6, 534}}};

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
        workerThread.lowPlyHistory[ss->ply][move.from_to()] << bonus * 879 / 1024;

    update_continuation_histories(ss, pos.moved_piece(move), move.to_sq(), bonus * 888 / 1024);

    int pIndex = pawn_structure_index(pos);
    workerThread.pawnHistory[pIndex][pos.moved_piece(move)][move.to_sq()] << bonus * 634 / 1024;
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
