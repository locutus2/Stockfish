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

#include "bitboard.h"
#include "evaluate.h"
#include "history.h"
#include "misc.h"
#include "movegen.h"
#include "movepick.h"
#include "nnue/network.h"
#include "nnue/nnue_accumulator.h"
#include "position.h"
#include "syzygy/tbprobe.h"
#include "thread.h"
#include "timeman.h"
#include "tt.h"
#include "types.h"
#include "uci.h"
#include "ucioption.h"

namespace Stockfish {

static constexpr std::array<int, 16> lmrDivisor = {3307, 2930, 2874, 2818, 3215, 3225, 3224, 2782,
                                                   2858, 2919, 3088, 3275, 3180, 2868, 3006, 3599};

namespace TB = Tablebases;

void syzygy_extend_pv(const OptionsMap&            options,
                      const Search::LimitsType&    limits,
                      Stockfish::Position&         pos,
                      Stockfish::Search::RootMove& rootMove,
                      Value&                       v);

using namespace Search;

namespace {

constexpr u64 NODES_LIMIT_OUTPUT = 10'000'000;

constexpr int SEARCHEDLIST_CAPACITY = 32;
using SearchedList                  = ValueList<Move, SEARCHEDLIST_CAPACITY>;

// (*Scalers):
// The values with Scaler asterisks have proven non-linear scaling.
// They are optimized to time controls of 180 + 1.8 and longer,
// so changing them or adding conditions that are similar requires
// tests at these types of time controls.

// (*Scaler) All tuned parameters at time controls shorter than
// optimized for require verifications at longer time controls

int correction_value(const Worker& w, const Position& pos, const Stack* const ss) {
    const Color us     = pos.side_to_move();
    const auto  m      = (ss - 1)->currentMove;
    const auto& shared = w.sharedHistory;
    const int   pcv    = shared.pawn_correction_entry(pos)[us].pawn;
    const int   micv   = shared.minor_piece_correction_entry(pos)[us].minor;
    const int   wnpcv  = shared.nonpawn_correction_entry<WHITE>(pos)[us].nonPawnWhite;
    const int   bnpcv  = shared.nonpawn_correction_entry<BLACK>(pos)[us].nonPawnBlack;
    const int   cntcv =
      m.is_ok()
          ? 8363
            * ((*(ss - 2)->continuationCorrectionHistory)[pos.piece_on(m.to_sq())][m.to_sq()]
               + (*(ss - 4)->continuationCorrectionHistory)[pos.piece_on(m.to_sq())][m.to_sq()])
          : 64549;

    return 13345 * pcv + 9280 * micv + 11840 * (wnpcv + bnpcv) + cntcv;
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

    constexpr int nonPawnWeight = 186;
    auto&         shared        = workerThread.sharedHistory;

    shared.pawn_correction_entry(pos)[us].pawn << bonus;
    shared.minor_piece_correction_entry(pos)[us].minor << bonus * 152 / 128;
    shared.nonpawn_correction_entry<WHITE>(pos)[us].nonPawnWhite << bonus * nonPawnWeight / 128;
    shared.nonpawn_correction_entry<BLACK>(pos)[us].nonPawnBlack << bonus * nonPawnWeight / 128;

    if (m.is_ok())
    {
        const Square to = m.to_sq();
        const Piece  pc = pos.piece_on(to);
        (*(ss - 2)->continuationCorrectionHistory)[pc][to] << bonus * 136 / 128;
        (*(ss - 4)->continuationCorrectionHistory)[pc][to] << bonus * 68 / 128;
    }
}

// Add a small random component to draw evaluations to avoid 3-fold blindness
Value value_draw(usize nodes) { return VALUE_DRAW - 1 + Value(nodes & 0x2); }
Value value_to_tt(Value v, int ply);
Value value_from_tt(Value v, int ply, int r50c);
void  update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus);
void  update_quiet_histories(
   const Position& pos, Stack* ss, Search::Worker& workerThread, Move move, int bonus);
void update_all_stats(const Position& pos,
                      Stack*          ss,
                      Search::Worker& workerThread,
                      Move            bestMove,
                      Square          prevSq,
                      SearchedList&   quietsSearched,
                      SearchedList&   capturesSearched,
                      Depth           depth,
                      Move            ttMove,
                      bool            PvNode);

// Detect shuffling moves in order to limit search explosions
// Added in #6447 as non-regression, and so its parameters should not be tuned
bool is_shuffling(Move move, Stack* const ss, const Position& pos) {
    if (pos.capture_stage(move) || pos.rule50_count() < 10)
        return false;
    if (pos.state()->pliesFromNull < 6 || ss->ply < 20)
        return false;
    return move.from_sq() == (ss - 2)->currentMove.to_sq()
        && (ss - 2)->currentMove.from_sq() == (ss - 4)->currentMove.to_sq();
}

}  // namespace

Search::Worker::Worker(SharedState&                    sharedState,
                       std::unique_ptr<ISearchManager> sm,
                       usize                           threadId,
                       usize                           numaThreadId,
                       usize                           numaTotalThreads,
                       NumaReplicatedAccessToken       token) :
    // Unpack the SharedState struct into member variables
    sharedHistory(sharedState.sharedHistories.at(token.get_numa_index())),
    continuationHistory(sharedHistory.continuationHistory),
    threadIdx(threadId),
    numaThreadIdx(numaThreadId),
    numaTotal(numaTotalThreads),
    numaAccessToken(token),
    manager(std::move(sm)),
    options(sharedState.options),
    threads(sharedState.threads),
    tt(sharedState.tt),
    network(sharedState.network),
    refreshTable(network[token]) {
    clear();
}

void Search::Worker::ensure_network_replicated() {
    // Access once to force lazy initialization.
    // We do this because we want to avoid initialization during search.
    (void) (network[numaAccessToken]);
}

void Search::Worker::start_searching() {

    accumulatorStack.reset();

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
        main_manager()->updates.onUpdateNoMoves(
          {0, {rootPos.checkers() ? -VALUE_MATE : VALUE_DRAW, rootPos}});
        main_manager()->updates.onBestmove(UCIEngine::move(Move::none()), "");
        return;
    }

    // Main thread starts non-main threads, and begins own search.
    threads.start_searching();
    bool uciPvSent = iterative_deepening();

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

    if (!limits.depth && !skill.enabled())
        bestThread = threads.get_best_thread()->worker.get();

    main_manager()->bestPreviousScore        = bestThread->rootMoves[0].score;
    main_manager()->bestPreviousAverageScore = bestThread->rootMoves[0].averageScore;

    if (bestThread->rootMoves[0].pv.size() == 1
        && bestThread->rootMoves[0].extract_ponder_from_tt(tt, rootPos))
        uciPvSent = false;

    // Send PV info if it has changed since last output in iterative_deepening().
    if (!uciPvSent || bestThread != this)
        main_manager()->output_pv(*bestThread, threads, tt, bestThread->rootDepth);

    // In rare cases, output_pv() may change the ponder move through syzygy_extend_pv().
    std::string ponder;
    if (bestThread->rootMoves[0].pv.size() > 1)
        ponder = UCIEngine::move(bestThread->rootMoves[0].pv[1], rootPos.is_chess960());

    auto bestmove = UCIEngine::move(bestThread->rootMoves[0].pv[0], rootPos.is_chess960());
    main_manager()->updates.onBestmove(bestmove, ponder);
}

// Main iterative deepening loop. It calls search()
// repeatedly with increasing depth until the allocated thinking time has been
// consumed, the user stops the search, or the maximum search depth is reached.
bool Search::Worker::iterative_deepening() {

    SearchManager* mainThread = (is_mainthread() ? main_manager() : nullptr);

    PVMoves pv;

    PVMoves lastBestMovePV;
    Depth   lastBestMoveDepth = 0;
    Value   lastBestMoveScore = -VALUE_INFINITE;

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
          &continuationHistory[0][0][NO_PIECE][0];  // Use as a sentinel
        (ss - i)->continuationCorrectionHistory = &continuationCorrectionHistory[NO_PIECE][0];
        (ss - i)->staticEval                    = VALUE_NONE;
    }

    for (int i = 0; i <= MAX_PLY + 2; ++i)
        (ss + i)->ply = i;

    ss->pv = &pv;

    if (mainThread)
    {
        if (mainThread->bestPreviousScore == VALUE_INFINITE)
            mainThread->iterValue.fill(VALUE_ZERO);
        else
            mainThread->iterValue.fill(mainThread->bestPreviousScore);
    }

    usize multiPV = usize(options["MultiPV"]);
    Skill skill(options["Skill Level"], options["UCI_LimitStrength"] ? int(options["UCI_Elo"]) : 0);

    // When playing with strength handicap enable MultiPV search that we will
    // use behind-the-scenes to retrieve a set of possible moves.
    if (skill.enabled())
        multiPV = std::max(multiPV, usize(4));

    multiPV = std::min(multiPV, rootMoves.size());

    int  searchAgainCounter = 0;
    bool uciPvSent          = false;

    lowPlyHistory.fill(100);

    for (Color c : {WHITE, BLACK})
        for (int i = 0; i < UINT_16_HISTORY_SIZE; i++)
            mainHistory[c][i] = mainHistory[c][i] * 789 / 1024;

    // Iterative deepening loop until requested to stop or the target depth is reached
    while (rootDepth + 1 < MAX_PLY && !threads.stop
           && !(limits.depth && mainThread && rootDepth >= limits.depth))
    {
        rootDepth++;

        // Age out PV variability metric and signal the start of a new iteration.
        if (mainThread)
        {
            totBestMoveChanges /= 2;
            uciPvSent = false;
        }

        // Save the last iteration's scores before the first PV line is searched and
        // all the move scores except the (new) PV are set to -VALUE_INFINITE.
        for (usize i = 0; i < rootMoves.size(); ++i)
        {
            rootMoves[i].previousScore      = rootMoves[i].score;
            rootMoves[i].previousPV         = rootMoves[i].pv;
            rootMoves[i].previousScoreExact = i < multiPV;
        }

        usize pvFirst = pvLast = 0;

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

            lastIterationIdxPV = rootMoves[pvIdx].previousPV;

            // Reset UCI info selDepth for each depth and each PV line
            selDepth = 0;

            // Reset aspiration window starting size
            delta     = 5 + threadIdx % 8 + std::abs(rootMoves[pvIdx].meanSquaredScore) / 10588;
            Value avg = rootMoves[pvIdx].averageScore;
            alpha     = std::max(avg - delta, -VALUE_INFINITE);
            beta      = std::min(avg + delta, VALUE_INFINITE);

            // Adjust optimism based on root move's averageScore
            optimism[us]  = 137 * avg / (std::abs(avg) + 81);
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
                    && nodes > NODES_LIMIT_OUTPUT)
                    main_manager()->output_pv(*this, threads, tt, rootDepth);

                // In case of failing low/high increase aspiration window and re-search,
                // otherwise exit the loop.
                if (bestValue <= alpha)
                {
                    beta  = alpha;
                    alpha = std::max(bestValue - delta, -VALUE_INFINITE);

                    failedHighCnt = 0;
                    if (mainThread)
                        mainThread->stopOnPonderhit = false;
                }
                else if (bestValue >= beta)
                {
                    alpha = std::max(beta - delta, alpha);
                    beta  = std::min(bestValue + delta, VALUE_INFINITE);
                    ++failedHighCnt;
                }
                else
                    break;

                delta += 44 * delta / 128;

                assert(alpha >= -VALUE_INFINITE && beta <= VALUE_INFINITE);
            }

            if (threads.stop && pvIdx)
            {
                // In multiPV analysis we do not let aborted searches spoil mated-in/
                // TB loss scores from a completed search in an earlier PV line.
                // Hence we guard against an aborted pvIdx line overtaking pvIdx - 1
                // when pvIdx - 1 is a proven loss.
                // Moreover, we do not trust an exact loss score from an aborted search.
                if ((is_loss(rootMoves[pvIdx - 1].score) && rootMoves[pvIdx] < rootMoves[pvIdx - 1])
                    || rootMoves[pvIdx].score_is_exact_loss())
                {
                    // If previousScore is exact and worse than pvIdx - 1, we can safely use it.
                    // If it is equal, we make sure it cannot overtake pvIdx - 1.
                    if (rootMoves[pvIdx].previousScore != -VALUE_INFINITE
                        && rootMoves[pvIdx].previousScoreExact
                        && rootMoves[pvIdx].previousScore <= rootMoves[pvIdx - 1].score)
                    {
                        rootMoves[pvIdx].score = rootMoves[pvIdx].uciScore =
                          rootMoves[pvIdx].previousScore;
                        rootMoves[pvIdx].previousScore = -VALUE_INFINITE;
                        rootMoves[pvIdx].pv            = rootMoves[pvIdx].previousPV;
                        rootMoves[pvIdx].unset_bound_flags();
                    }

                    // Otherwise, if we can, we cap the score to the best possible, and mark
                    // the score as a bound (also a valid excuse for the incomplete PV.)
                    else
                    {
                        if (is_loss(rootMoves[pvIdx - 1].score))
                        {
                            rootMoves[pvIdx].score = rootMoves[pvIdx].uciScore =
                              rootMoves[pvIdx - 1].score;
                            rootMoves[pvIdx].previousScore = -VALUE_INFINITE;
                            rootMoves[pvIdx].pv.resize(1);
                            rootMoves[pvIdx].scoreUpperbound = true;
                        }
                        else
                            rootMoves[pvIdx].scoreUpperbound = false;

                        rootMoves[pvIdx].scoreLowerbound = !rootMoves[pvIdx].scoreUpperbound;
                    }
                }

                // Finally, we mark all loss scores from partially searched moves as a bound.
                for (usize i = pvIdx + 1; i < multiPV; ++i)
                    if (rootMoves[i].score_is_exact_loss())
                        rootMoves[i].scoreLowerbound = true;
            }

            // Sort the PV lines searched so far and update the GUI
            std::stable_sort(rootMoves.begin() + pvFirst, rootMoves.begin() + pvIdx + 1);

            if (mainThread && !threads.stop && (pvIdx + 1 == multiPV || nodes > NODES_LIMIT_OUTPUT))
            {
                main_manager()->output_pv(*this, threads, tt, rootDepth);
                uciPvSent = (pvIdx + 1 == multiPV);
            }

            if (threads.stop)
                break;
        }

        const bool forgottenMate = lastBestMoveScore != -VALUE_INFINITE
                                && is_mate_or_mated(lastBestMoveScore)
                                && (std::abs(rootMoves[0].score) < std::abs(lastBestMoveScore)
                                    || rootMoves[0].score_is_bound());

        if (!threads.stop)
        {
            if (lastBestMovePV.empty() || lastBestMovePV[0] != rootMoves[0].pv[0])
                lastBestMoveDepth = rootDepth;

            // Do not replace (shorter) mate scores from a previous iteration.
            if (!forgottenMate)
            {
                lastBestMovePV    = rootMoves[0].pv;
                lastBestMoveScore = rootMoves[0].score;
            }
        }

        const bool abortedLossSearch = threads.stop && !pvIdx && rootMoves[0].score_is_exact_loss();

        // An exact mated-in/TB-loss score from an aborted search cannot be trusted: the
        // loss could be delayed or refuted upon exploring the remaining root-moves.
        // Thus here we roll back to the score from the previous iteration.
        // We do the same if a search has failed to recover a mate score that was found
        // in a previous iteration.
        if (abortedLossSearch || (rootMoves[0].score != -VALUE_INFINITE && forgottenMate))
        {
            // Bring the last best move to the front for best thread selection.
            if (!lastBestMovePV.empty())
            {
                Utility::move_to_front(rootMoves, [&lastPV = std::as_const(lastBestMovePV)](
                                                    const auto& rm) { return rm == lastPV[0]; });
                rootMoves[0].score = rootMoves[0].uciScore = lastBestMoveScore;
                rootMoves[0].pv                            = lastBestMovePV;
                rootMoves[0].unset_bound_flags();

                if (mainThread)
                    uciPvSent = false;
            }
            // For an aborted d1 search we label the loss score as a lower bound.
            else if (abortedLossSearch)
                rootMoves[0].scoreLowerbound = true;
        }

        // Have we found a "mate in x" after a completed iteration?
        if (limits.mate && !threads.stop && is_mate_or_mated(rootMoves[0].score)
            && VALUE_MATE - std::abs(rootMoves[0].score) <= 2 * limits.mate)
            threads.stop = true;

        if (!mainThread)
            continue;

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
            u64 nodesEffort = rootMoves[0].effort * 100000 / std::max(u64(1), u64(nodes));

            double fallingEval = (11.87 + 2.21 * (mainThread->bestPreviousAverageScore - bestValue)
                                  + 1.0 * (mainThread->iterValue[iterIdx] - bestValue))
                               / 100.0;
            fallingEval = std::clamp(fallingEval, 0.572, 1.708);

            // If the bestMove is stable over several iterations, reduce time accordingly
            timeReduction =
              std::clamp(interpolate(double(rootDepth - lastBestMoveDepth), 5.0, 18.0, 0.65, 1.55),
                         0.65, 1.55);

            double reduction = (1.48 + mainThread->previousTimeReduction) / (2.157 * timeReduction);

            double bestMoveInstability = 1.096 + 2.29 * totBestMoveChanges / threads.size();

            double highBestMoveEffort = std::clamp(
              interpolate(i64(nodesEffort), i64(79219), i64(101822), 0.924, 0.71), 0.71, 0.924);

            double totalTime = mainThread->tm.optimum() * fallingEval * reduction
                             * bestMoveInstability * highBestMoveEffort;

            if (rootMoves.size() == 1)
                // Cap used time to 0.5s for a better viewer experience
                totalTime = std::min(500.0, totalTime);

            auto elapsedTime = elapsed();

            // Stop the search if we have exceeded totalTime or maximum time,
            // or if we know that there are no better moves in the analysed line(s)
            if (elapsedTime > std::min(totalTime, double(mainThread->tm.maximum()))
                || rootMoves[multiPV - 1].score >= mate_in(3) || rootMoves[0].score == mated_in(2))
            {
                // If we are allowed to ponder do not stop the search now but
                // keep pondering until the GUI sends "ponderhit" or "stop".
                if (mainThread->ponder)
                    mainThread->stopOnPonderhit = true;
                else
                    threads.stop = true;
            }
            else
                threads.increaseDepth = mainThread->ponder || elapsedTime <= totalTime * 0.50;
        }

        mainThread->iterValue[iterIdx] = bestValue;
        iterIdx                        = (iterIdx + 1) & 3;
    }

    if (!mainThread)
        return false;

    mainThread->previousTimeReduction = timeReduction;

    // If the skill level is enabled, swap the best PV line with the sub-optimal one
    if (skill.enabled())
        std::swap(rootMoves[0],
                  *std::find(rootMoves.begin(), rootMoves.end(),
                             skill.best ? skill.best : skill.pick_best(rootMoves, multiPV)));

    return uciPvSent;
}


void Search::Worker::do_move(Position& pos, const Move move, StateInfo& st, Stack* const ss) {
    do_move(pos, move, st, pos.gives_check(move), ss);
}

void Search::Worker::do_move(
  Position& pos, const Move move, StateInfo& st, const bool givesCheck, Stack* const ss) {
    // prefetch_key does not model castling, en passant or promotion keys
    // exactly; for rare moves the prefetch lands on an unused line.
    prefetch(tt.first_entry(pos.prefetch_key(move)));

    bool capture = pos.capture_stage(move);
    ++nodes;

    auto [dirtyPiece, dirtyThreats] = accumulatorStack.push();
    pos.do_move(move, st, givesCheck, dirtyPiece, dirtyThreats, &tt, &sharedHistory);

    if (ss != nullptr)
    {
        ss->currentMove = move;
        ss->continuationHistory =
          &continuationHistory[ss->inCheck][capture][dirtyPiece.pc][move.to_sq()];
        ss->continuationCorrectionHistory =
          &continuationCorrectionHistory[dirtyPiece.pc][move.to_sq()];
    }
}

void Search::Worker::do_null_move(Position& pos, StateInfo& st, Stack* const ss) {
    pos.do_null_move(st);
    ss->currentMove                   = Move::null();
    ss->continuationHistory           = &continuationHistory[0][0][NO_PIECE][0];
    ss->continuationCorrectionHistory = &continuationCorrectionHistory[NO_PIECE][0];
}

void Search::Worker::undo_move(Position& pos, const Move move) {
    pos.undo_move(move);
    accumulatorStack.pop();
}

void Search::Worker::undo_null_move(Position& pos) { pos.undo_null_move(); }


// Reset histories, usually before a new game
void Search::Worker::clear() {
    mainHistory.fill(-5);
    captureHistory.fill(-699);

    // Each thread is responsible for clearing their part of shared history
    sharedHistory.correctionHistory.clear_range(-6, numaThreadIdx, numaTotal);
    sharedHistory.pawnHistory.clear_range(-1262, numaThreadIdx, numaTotal);

    ttMoveHistory = 0;

    for (auto& to : continuationCorrectionHistory)
        for (auto& h : to)
            h.fill(5);

    for (bool inCheck : {false, true})
        for (StatsType c : {NoCaptures, Captures})
            for (auto& to : continuationHistory[inCheck][c])
                for (auto& h : to)
                    h.fill(-552);

    for (usize i = 1; i < reductions.size(); ++i)
        reductions[i] = int(2834 / 128.0 * std::log(i));

    refreshTable.clear(network[numaAccessToken]);
}


// Main search function for both PV and non-PV nodes
template<NodeType nodeType>
Value Search::Worker::search(
  Position& pos, Stack* ss, Value alpha, Value beta, Depth depth, const bool cutNode) {

    constexpr bool PvNode   = nodeType != NonPV;
    constexpr bool rootNode = nodeType == Root;
    const bool     allNode  = !(PvNode || cutNode);

    // Dive into quiescence search when the depth reaches zero
    if (depth <= 0)
        return qsearch<PvNode ? PV : NonPV>(pos, ss, alpha, beta);

    bool CC = false;
    bool C = true;
    Value V = 0;
    Depth origDepth = depth;
    Value origAlpha = alpha;
    if(depth < 32)
    {
	    for(int d = 0; d < 32; d++)
		    dbg_hit_on(d == depth, d);
    }

    // Limit the depth if extensions made it too large
    depth = std::min(depth, MAX_PLY - 1);

    // Check if we have an upcoming move that draws by repetition
    if (!rootNode && alpha < VALUE_DRAW && pos.upcoming_repetition(ss->ply))
    {
        alpha = value_draw(nodes);
        if (alpha >= beta)
	{
	    if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	    if(origDepth < 32) dbg_hit_on(false, 200 + origDepth);
            return alpha;
	}
    }

    assert(-VALUE_INFINITE <= alpha && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - 1));
    assert(0 < depth && depth < MAX_PLY);
    assert(!(PvNode && cutNode));

    PVMoves   pv;
    StateInfo st;

    Key   posKey;
    Move  move, excludedMove, bestMove;
    Depth extension, newDepth;
    Value bestValue, value, eval, maxValue, probCutBeta;
    bool  givesCheck, improving, priorCapture, opponentWorsening;
    bool  capture, ttCapture;
    int   priorReduction;
    Piece movedPiece;

    SearchedList capturesSearched;
    SearchedList quietsSearched;

    // Step 1. Initialize node
    ss->inCheck   = pos.checkers();
    priorCapture  = pos.captured_piece();
    Color us      = pos.side_to_move();
    ss->moveCount = 0;
    bestValue     = -VALUE_INFINITE;
    maxValue      = VALUE_INFINITE;

    ss->followPV = rootNode
                || ((ss - 1)->followPV
                    && (static_cast<usize>(ss->ply - 1) < lastIterationIdxPV.size()
                        && (ss - 1)->currentMove == lastIterationIdxPV[ss->ply - 1]));

    // Check for the available remaining time
    if (is_mainthread())
        main_manager()->check_time(*this);

    // Used to send selDepth info to GUI (selDepth counts from 1, ply from 0)
    if (PvNode && selDepth < ss->ply + 1)
        selDepth = ss->ply + 1;

    if (!rootNode)
    {
        // Step 2. Check for aborted search and immediate draw
        if (threads.stop.load(std::memory_order_relaxed) || pos.is_draw(ss->ply)
            || ss->ply >= MAX_PLY)
	{
	    if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	    if(origDepth < 32) dbg_hit_on(((ss->ply >= MAX_PLY && !ss->inCheck) ? evaluate(pos) : value_draw(nodes)) > alpha, 200 + origDepth);
            return (ss->ply >= MAX_PLY && !ss->inCheck) ? evaluate(pos) : value_draw(nodes);
	}

        // Step 3. Mate distance pruning. Even if we mate at the next move our score
        // would be at best mate_in(ss->ply + 1), but if alpha is already bigger because
        // a shorter mate was found upward in the tree then there is no need to search
        // because we will never beat the current alpha. Same logic but with reversed
        // signs apply also in the opposite condition of being mated instead of giving
        // mate. In this case, return a fail-high score.
        alpha = std::max(mated_in(ss->ply), alpha);
        beta  = std::min(mate_in(ss->ply + 1), beta);
        if (alpha >= beta)
	{
	    if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	    if(origDepth < 32) dbg_hit_on(false, 200 + origDepth);
            return alpha;
	}
    }

    assert(0 <= ss->ply && ss->ply < MAX_PLY);

    Square prevSq  = ((ss - 1)->currentMove).is_ok() ? ((ss - 1)->currentMove).to_sq() : SQ_NONE;
    bestMove       = Move::none();
    priorReduction = (ss - 1)->reduction;
    (ss - 1)->reduction = 0;
    ss->statScore       = 0;
    (ss + 2)->cutoffCnt = 0;

    const auto correctionValue = correction_value(*this, pos, ss);

    // Step 4. Transposition table lookup
    excludedMove                   = ss->excludedMove;
    posKey                         = pos.key();
    auto [ttHit, ttData, ttWriter] = tt.probe(posKey);
    // Need further processing of the saved data
    ss->ttHit    = ttHit;
    ttData.move  = rootNode ? rootMoves[pvIdx].pv[0] : ttHit ? ttData.move : Move::none();
    ttData.value = ttHit ? value_from_tt(ttData.value, ss->ply, pos.rule50_count()) : VALUE_NONE;
    ss->ttPv     = excludedMove ? ss->ttPv : PvNode || (ttHit && ttData.is_pv);
    ttCapture    = ttData.move && pos.capture_stage(ttData.move);

    // Step 5. Static evaluation of the position
    Value unadjustedStaticEval = VALUE_NONE;

    // Skip early pruning when in check
    if (ss->inCheck)
        ss->staticEval = eval = (ss - 2)->staticEval;
    else if (excludedMove)
        unadjustedStaticEval = eval = ss->staticEval;
    else if (ss->ttHit)
    {
        // Never assume anything about values stored in TT
        unadjustedStaticEval = ttData.eval;
        if (!is_valid(unadjustedStaticEval))
            unadjustedStaticEval = evaluate(pos);

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

    // Set up the improving flag, which is true if current static evaluation is
    // bigger than the previous static evaluation at our turn (if we were in
    // check at our previous move we go back until we weren't in check) and is
    // false otherwise. The improving flag is used in various pruning heuristics.
    // Similarly, opponentWorsening is true if our static evaluation is better
    // for us than at the last ply.
    improving         = ss->staticEval > (ss - 2)->staticEval;
    opponentWorsening = ss->staticEval > -(ss - 1)->staticEval;

    // Hindsight adjustment of reductions based on static evaluation difference.
    if (priorReduction >= 3 && !opponentWorsening)
        depth++;
    if (priorReduction >= 2 && depth >= 2 && ss->staticEval + (ss - 1)->staticEval > 173)
        depth--;

    // At non-PV nodes we check for an early TT cutoff
    if (!PvNode && !excludedMove && ttData.depth > depth - (ttData.value <= beta)
        && is_valid(ttData.value)  // Can happen when !ttHit or when access race in probe()
        && (ttData.bound & (ttData.value >= beta ? BOUND_LOWER : BOUND_UPPER))
        && (cutNode == (ttData.value >= beta) || depth > 4))
    {
        // If ttMove is quiet, update move sorting heuristics on TT hit
        if (ttData.move && ttData.value >= beta)
        {
            // Bonus for a quiet ttMove that fails high
            if (!ttCapture)
                update_quiet_histories(pos, ss, *this, ttData.move, std::min(114 * depth, 724));

            // Extra penalty for early quiet moves of the previous ply
            if (prevSq != SQ_NONE && (ss - 1)->moveCount < 4 && !priorCapture)
                update_continuation_histories(ss - 1, pos.piece_on(prevSq), prevSq, -2187);
        }

        // Partial workaround for the graph history interaction problem
        // For high rule50 counts don't produce transposition table cutoffs.
        if (pos.rule50_count() < 96)
        {
            if (depth >= 7 && ttData.move && pos.pseudo_legal(ttData.move) && pos.legal(ttData.move)
                && !is_decisive(ttData.value))
            {
                pos.do_move(ttData.move, st);
                Key nextPosKey                             = pos.key();
                auto [ttHitNext, ttDataNext, ttWriterNext] = tt.probe(nextPosKey);
                pos.undo_move(ttData.move);

                // Check that the ttValue after the tt move would also trigger a cutoff
                if (!is_valid(ttDataNext.value))
		{
	            if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	            if(origDepth < 32) dbg_hit_on(ttData.value > alpha, 200 + origDepth);
                    return ttData.value;
		}

                if ((ttData.value >= beta) == (-ttDataNext.value >= beta))
		{
	            if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	            if(origDepth < 32) dbg_hit_on(ttData.value > alpha, 200 + origDepth);
                    return ttData.value;
		}
            }
            else
	    {
	        if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	        if(origDepth < 32) dbg_hit_on(ttData.value > alpha, 200 + origDepth);
                return ttData.value;
	    }
        }
    }  // No cutoff, but why? Does the stored inexact value mismatch our aspiration window?
    else if (!PvNode && !excludedMove && ttData.depth > depth - (ttData.value <= beta)
             && is_valid(ttData.value) && ttData.bound != BOUND_EXACT
             && ttData.bound & (ttData.value >= beta ? BOUND_UPPER : BOUND_LOWER) && depth > 5)
    {  // If a window-bound mismatch is the only reason cutoff failed, penalize the now-useless tte
        ttWriter.penalize(1);
    }

    // Step 6. Tablebases probe
    if (!rootNode && !excludedMove && tbConfig.cardinality)
    {
        int piecesCount = pos.count<ALL_PIECES>();

        if (piecesCount <= tbConfig.cardinality
            && (piecesCount < tbConfig.cardinality || depth >= tbConfig.probeDepth)
            && pos.rule50_count() == 0 && !pos.can_castle(ANY_CASTLING))
        {
            TB::ProbeState err;
            TB::WDLScore   wdl = TB::probe_wdl(pos, &err);

            // Force check of time on the next occasion
            if (is_mainthread())
                main_manager()->callsCnt = 0;

            if (err != TB::ProbeState::FAIL)
            {
                ++tbHits;

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

    if (ss->inCheck)
        goto moves_loop;

    // Use static evaluation difference to improve quiet move ordering
    if (((ss - 1)->currentMove).is_ok() && !(ss - 1)->inCheck && !priorCapture)
    {
        int evalDiff = std::clamp(-int((ss - 1)->staticEval + ss->staticEval), -183, 180) + 62;
        mainHistory[~us][((ss - 1)->currentMove).raw()] << evalDiff * 10;
        if (!ttHit && type_of(pos.piece_on(prevSq)) != PAWN
            && ((ss - 1)->currentMove).type_of() != PROMOTION)
            sharedHistory.pawn_entry(pos)[pos.piece_on(prevSq)][prevSq] << evalDiff * 13;
    }


    // Step 7. Razoring
    // If eval is really low, skip search entirely and return the qsearch value.
    // For PvNodes, we must have a guard against mates being returned.
    if (!PvNode && eval < alpha - 465 - 300 * depth * depth)
    {
	CC = true;
	C = bool(excludedMove);
	if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
        V = qsearch<NonPV>(pos, ss, alpha, beta);
	if(origDepth < 32) dbg_hit_on(V > alpha, 200 + origDepth);
        if(!CC) return V;
    }

    // Step 8. Futility pruning: child node
    // The depth condition is important for mate finding.
    if (!ss->ttPv && depth < 17 && eval >= beta && (!ttData.move || ttCapture) && !is_loss(beta)
        && !is_win(eval))
    {
        Value futilityMult = std::min(40 + depth * 4, 80);
        futilityMult -= 20 * !ss->ttHit;

        Value futilityMargin = futilityMult * depth
                             - (2934 * improving + 343 * opponentWorsening) * futilityMult / 1024
                             + std::abs(correctionValue) / 182069;

        if (eval - futilityMargin >= beta)
	{
	    if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	    if(origDepth < 32) dbg_hit_on((716 * beta + 308 * eval) / 1024 > alpha, 200 + origDepth);
            return (716 * beta + 308 * eval) / 1024;
	}
    }

    // Step 9. Null move search with verification search
    if (cutNode && ss->staticEval >= beta - 14 * depth - 45 * improving + 374 && !excludedMove
        && pos.non_pawn_material(us) && ss->ply >= nmpMinPly && !is_loss(beta))
    {
        assert((ss - 1)->currentMove != Move::null());

        // Null move dynamic reduction based on depth
        Depth R = 7 + depth / 3;
        do_null_move(pos, st, ss);

        Value nullValue = -search<NonPV>(pos, ss + 1, -beta, -beta + 1, depth - R, false);

        undo_null_move(pos);

        // Do not return unproven mate or TB scores
        if (nullValue >= beta && !is_win(nullValue))
        {
            if (nmpMinPly || depth < 16)
	    {
	        if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	        if(origDepth < 32) dbg_hit_on(nullValue > alpha, 200 + origDepth);
                return nullValue;
	    }

            assert(!nmpMinPly);  // Recursive verification is not allowed

            // Do verification search at high depths, with null move pruning disabled
            // until ply exceeds nmpMinPly.
            nmpMinPly = ss->ply + 3 * (depth - R) / 4;

            Value v = search<NonPV>(pos, ss, beta - 1, beta, depth - R, false);

            nmpMinPly = 0;

            if (v >= beta)
	    {
	        if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	        if(origDepth < 32) dbg_hit_on(nullValue > alpha, 200 + origDepth);
                return nullValue;
	    }
        }
    }

    improving |= ss->staticEval >= beta;

    // Step 10. Internal iterative reductions
    // At sufficient depth, reduce depth for PV/Cut nodes without a TTMove.
    // (*Scaler) Making IIR more aggressive scales poorly.
    if (!ss->followPV && !allNode && depth >= 6 && !ttData.move)
        depth--;

    // Step 11. ProbCut
    // If we have a good enough capture (or queen promotion) and a reduced search
    // returns a value much above beta, we can (almost) safely prune the previous move.
    probCutBeta = beta + 214 - 59 * improving;
    if (depth >= 3
        && !is_decisive(beta)
        // If value from transposition table is lower than probCutBeta, don't attempt
        // probCut there
        && !(is_valid(ttData.value) && ttData.value < probCutBeta))
    {
        assert(probCutBeta < VALUE_INFINITE && probCutBeta > beta);

        MovePicker mp(pos, ttData.move, probCutBeta - ss->staticEval, &captureHistory);
        Depth      probCutDepth = depth - 4 - improving;

        while ((move = mp.next_move()) != Move::none())
        {
            assert(move.is_ok());

            if (move == excludedMove || !pos.legal(move))
                continue;

            assert(pos.capture_stage(move));

            do_move(pos, move, st, ss);

            // Perform a preliminary qsearch to verify that the move holds
            value = -qsearch<NonPV>(pos, ss + 1, -probCutBeta, -probCutBeta + 1);

            // If the qsearch held, perform the regular search
            if (value >= probCutBeta && probCutDepth > 0)
                value = -search<NonPV>(pos, ss + 1, -probCutBeta, -probCutBeta + 1, probCutDepth,
                                       !cutNode);

            undo_move(pos, move);

            if (value >= probCutBeta)
            {
                // Save ProbCut data into transposition table
                ttWriter.write(posKey, value_to_tt(value, ss->ply), ss->ttPv, BOUND_LOWER,
                               probCutDepth + 1, move, unadjustedStaticEval, tt.generation());

                if (!is_decisive(value))
		{
	            if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	            if(origDepth < 32) dbg_hit_on(value - (probCutBeta - beta) > alpha, 200 + origDepth);
                    return value - (probCutBeta - beta);
		}
            }
        }
    }

moves_loop:  // When in check, search starts here

    // Step 12. A small Probcut idea
    probCutBeta = beta + 428;
    if ((ttData.bound & BOUND_LOWER) && ttData.depth >= depth - 4 && ttData.value >= probCutBeta
        && !is_decisive(beta) && is_valid(ttData.value) && !is_decisive(ttData.value))
    {
	if(origDepth < 32) dbg_hit_on(true, 100 + origDepth);
	if(origDepth < 32) dbg_hit_on(probCutBeta > alpha, 200 + origDepth);
        return probCutBeta;
    }

    if(origDepth < 32) dbg_hit_on(false, 100 + origDepth);

    const PieceToHistory* contHist[] = {
      (ss - 1)->continuationHistory, (ss - 2)->continuationHistory, (ss - 3)->continuationHistory,
      (ss - 4)->continuationHistory, (ss - 5)->continuationHistory, (ss - 6)->continuationHistory};


    MovePicker mp(pos, ttData.move, depth, &mainHistory, &lowPlyHistory, &captureHistory, contHist,
                  &sharedHistory, ss->ply);

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
        if (rootNode && !std::count(rootMoves.begin() + pvIdx, rootMoves.begin() + pvLast, move))
            continue;

        ss->moveCount = ++moveCount;

        if (rootNode && is_mainthread() && nodes > NODES_LIMIT_OUTPUT)
        {
            main_manager()->updates.onIter(
              {depth, UCIEngine::move(move, pos.is_chess960()), moveCount + pvIdx});
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

        int r = reduction(improving, depth, moveCount, delta);

        // Increase reduction for ttPv nodes (*Scaler)
        // Larger values scale well
        if (ss->ttPv)
            r += 1006;

        // Step 14. Pruning at shallow depths.
        // Depth conditions are important for mate finding.
        if (!rootNode && pos.non_pawn_material(us) && !is_loss(bestValue))
        {
            // Skip quiet moves if movecount exceeds our threshold
            if (moveCount >= (3 + depth * depth) / (2 - improving))
                mp.skip_quiet_moves();

            // Reduced depth of the next LMR search
            int lmrDepth = newDepth - r / 1024;

            if (capture || givesCheck)
            {
                Piece capturedPiece = pos.piece_on(move.to_sq());
                int   captHist = captureHistory[movedPiece][move.to_sq()][type_of(capturedPiece)];

                // Futility pruning for captures
                if (!givesCheck && lmrDepth < 7)
                {
                    Value futilityValue = ss->staticEval + 231 + 232 * lmrDepth
                                        + PieceValue[capturedPiece] + 131 * captHist / 1024;

                    if (futilityValue <= alpha)
                        continue;
                }

                // SEE based pruning for captures and checks
                // Avoid pruning sacrifices of our last piece for stalemate
                int margin = 175 * depth + captHist * 34 / 1024;
                if ((alpha >= VALUE_DRAW || pos.non_pawn_material(us) != PieceValue[movedPiece])
                    && !pos.see_ge(move, -margin))
                    continue;
            }
            else if (!ss->followPV || !PvNode)
            {
                int dIndex  = std::min(int(depth), int(lmrDivisor.size())) - 1;
                int history = (*contHist[0])[movedPiece][move.to_sq()]
                            + (*contHist[1])[movedPiece][move.to_sq()]
                            + sharedHistory.pawn_entry(pos)[movedPiece][move.to_sq()];

                // Continuation history based pruning
                if (history < -4313 * depth)
                    continue;

                history += 64 * mainHistory[us][move.raw()] / 32;

                // (*Scaler): Generally, lower divisors scale well
                lmrDepth += history / lmrDivisor[dIndex];

                Value futilityValue = ss->staticEval + 40 + 138 * !bestMove + 117 * lmrDepth
                                    + 90 * (ss->staticEval > alpha);

                // Futility pruning: parent node
                // (*Scaler): Generally, more frequent futility pruning
                // scales well
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
        // Singular extension search. If all moves but one
        // fail low on a search of (alpha-s, beta-s), and just one fails high on
        // (alpha, beta), then that move is singular and should be extended. To
        // verify this we do a reduced search on the position excluding the ttMove
        // and if the result is lower than ttValue minus a margin, then we will
        // extend the ttMove. Recursive singular search is avoided.

        // (*Scaler) Generally, higher singularBeta (i.e closer to ttValue)
        // and lower extension margins scale well.
        if (!rootNode && move == ttData.move && !excludedMove && depth >= 6 + ss->ttPv
            && is_valid(ttData.value) && !is_decisive(ttData.value) && (ttData.bound & BOUND_LOWER)
            && ttData.depth >= depth - 3 && !is_shuffling(move, ss, pos))
        {
            Value singularBeta  = ttData.value - (60 + 70 * (ss->ttPv && !PvNode)) * depth / 59;
            Depth singularDepth = newDepth / 2;

            ss->excludedMove = move;
            value = search<NonPV>(pos, ss, singularBeta - 1, singularBeta, singularDepth, cutNode);
            ss->excludedMove = Move::none();

            if (value < singularBeta)
            {
                int corrValAdj   = std::abs(correctionValue) / 194822;
                int doubleMargin = -3 + 201 * PvNode - 157 * !ttCapture - corrValAdj
                                 - 1081 * ttMoveHistory / 117824 - (ss->ply > rootDepth) * 41;
                int tripleMargin = 72 + 306 * PvNode - 188 * !ttCapture + 84 * ss->ttPv - corrValAdj
                                 - (ss->ply > rootDepth) * 45;

                extension =
                  1 + (value < singularBeta - doubleMargin) + (value < singularBeta - tripleMargin);

                depth++;
            }

            // Multi-cut pruning
            // Our ttMove is assumed to fail high based on the bound of the TT entry,
            // and if after excluding the ttMove with a reduced search we fail high
            // over the original beta, we assume this expected cut-node is not
            // singular (multiple moves fail high), and we can prune the whole
            // subtree by returning a softbound.
            else if (value >= beta && !is_decisive(value))
            {
                ttMoveHistory << -442 - 108 * depth;
                return value;
            }

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

        u64 nodeCount = rootNode ? u64(nodes) : 0;

        // Step 16. Make the move
        do_move(pos, move, st, givesCheck, ss);

        // Add extension to new depth
        newDepth += extension;

        // Decrease reduction for PvNodes (*Scaler)
        if (ss->ttPv)
            r -= 2766 + PvNode * 1017 + (ttData.value > alpha) * 838
               + (ttData.depth >= depth) * (923 + cutNode * 955);

        r += 714;  // Base reduction offset to compensate for other tweaks
        r -= moveCount * 62;
        r -= std::abs(correctionValue) / 26131;

        // Increase reduction for cut nodes
        if (cutNode)
            r += 3995 + 1059 * !ttData.move;

        // Increase reduction if ttMove is a capture
        if (ttCapture)
            r += 1039;

        // Increase reduction if next ply has a lot of fail high
        if ((ss + 1)->cutoffCnt > 1)
            r += 236 + 1079 * ((ss + 1)->cutoffCnt > 2) + 1143 * allNode;

        // For first picked move (ttMove) reduce reduction
        else if (move == ttData.move)
            r -= 2016;

        if (capture)
            ss->statScore = 809 * int(PieceValue[pos.captured_piece()]) / 128
                          + captureHistory[movedPiece][move.to_sq()][type_of(pos.captured_piece())];
        else
            ss->statScore = 2 * mainHistory[us][move.raw()]
                          + (*contHist[0])[movedPiece][move.to_sq()]
                          + (*contHist[1])[movedPiece][move.to_sq()];

        // Decrease/increase reduction for moves with a good/bad history
        r -= ss->statScore * 445 / 4096;

        // Scale up reductions for expected ALL nodes
        if (allNode)
            r += r * 272 / (256 * depth + 285);

        // Step 17. Late moves reduction / extension (LMR)
        if (depth >= 2 && moveCount > 1)
        {
            // In general we want to cap the LMR depth search at newDepth, but when
            // reduction is negative, we allow this move a limited search extension
            // beyond the first move depth.
            // To prevent problems when the max value is less than the min value,
            // std::clamp has been replaced by a more robust implementation.
            Depth d = std::max(1, std::min(newDepth - r / 1024, newDepth + 2)) + PvNode;

            ss->reduction = newDepth - d;
            value         = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, d, true);
            ss->reduction = 0;

            // Do a full-depth search when reduced LMR search fails high
            // (*Scaler) Shallower searches here don't scale well
            if (value > alpha)
            {
                // Adjust full-depth search based on LMR results - if the result was
                // good enough search deeper, if it was bad enough search shallower.
                const bool doDeeperSearch    = d < newDepth && value > bestValue + 52;
                const bool doShallowerSearch = value < bestValue + 9;

                newDepth += doDeeperSearch - doShallowerSearch;

                if (newDepth > d)
                    value = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha, newDepth, !cutNode);

                // Post LMR continuation history updates
                update_continuation_histories(ss, movedPiece, move.to_sq(), 1415);
            }
        }

        // Step 18. Full-depth search when LMR is skipped
        else if (!PvNode || moveCount > 1)
        {
            // Increase reduction if ttMove is not present
            if (!ttData.move)
                r += 1085;

            // Note that if expected reduction is high, we reduce search depth here
            value = -search<NonPV>(pos, ss + 1, -(alpha + 1), -alpha,
                                   newDepth - (r > 5039) - (r > 5223 && newDepth > 2), !cutNode);
        }

        // For PV nodes only, do a full PV search on the first move or after a fail high,
        // otherwise let the parent node fail low with value <= alpha and try another move.
        if (PvNode && (moveCount == 1 || value > alpha))
        {
            (ss + 1)->pv = &pv;
            (ss + 1)->pv->clear();

            // Extend move from transposition table if we are about to dive into qsearch.
            // decisive score handling improves mate finding and retrograde analysis.
            if (move == ttData.move
                && ((is_valid(ttData.value) && is_decisive(ttData.value) && ttData.depth > 0)
                    || ttData.depth > 1))
                newDepth = std::max(newDepth, 1);

            value = -search<PV>(pos, ss + 1, -beta, -alpha, newDepth, false);
        }

        // Step 19. Undo move
        undo_move(pos, move);

        assert(value > -VALUE_INFINITE && value < VALUE_INFINITE);

        // Step 20. Check for a new best move
        // Finished searching the move. If a stop occurred, the return value of
        // the search cannot be trusted, and we return immediately without updating
        // best move, principal variation nor transposition table.
        if (threads.stop.load(std::memory_order_relaxed))
            return VALUE_ZERO;

        if (rootNode)
        {
            RootMove& rm = *std::find(rootMoves.begin(), rootMoves.end(), move);

            rm.effort += nodes - nodeCount;

            u64 N      = nodes - nodeCount;
            u64 E_prev = std::max(u64(1), rm.effort - N);

            // Dynamic EMA parameters for root move
            constexpr u64 Scale          = 32;
            constexpr u64 ChiNumerator   = 3;
            constexpr u64 ChiDenominator = 2;   // Chi = 3/2 = 1.5
            constexpr u64 MinWeight      = 12;  // 37.5% minimum weight
            constexpr u64 MaxWeight      = 24;  // 75% maximum weight

            u64 w     = std::clamp((Scale * N * ChiDenominator)
                                     / (N * ChiDenominator + ChiNumerator * E_prev),
                                   MinWeight, MaxWeight);
            u64 w_mss = std::min(w, u64(16));
            i64 v2    = i64(value) * std::abs(value);

            if (rm.averageScore == -VALUE_INFINITE)
                rm.averageScore = value;
            else
                rm.averageScore = Value((value * w + rm.averageScore * (Scale - w)) / Scale);

            if (rm.meanSquaredScore == -VALUE_INFINITE * VALUE_INFINITE)
                rm.meanSquaredScore = value * std::abs(value);
            else
                rm.meanSquaredScore =
                  Value((v2 * w_mss + int64_t(rm.meanSquaredScore) * (Scale - w_mss)) / Scale);

            // PV move or new best move?
            if (moveCount == 1 || value > alpha)
            {
                rm.score = rm.uciScore = value;
                rm.selDepth            = selDepth;
                rm.unset_bound_flags();

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

                for (Move pvMove : *(ss + 1)->pv)
                    rm.pv.push_back(pvMove);

                // We record how often the best move has been changed in each iteration.
                // This information is used for time management. In MultiPV mode,
                // we must take care to only do this for the first PV line.
                if (moveCount > 1 && !pvIdx)
                    ++bestMoveChanges;
            }
            else
                // All other moves but the PV, are set to the lowest value: this
                // is not a problem when sorting because the sort is stable and the
                // move position in the list is preserved - just the PV is pushed up.
                rm.score = -VALUE_INFINITE;
        }

        // In case we have an alternative move equal in eval to the current bestmove,
        // promote it to bestmove by pretending it just exceeds alpha (but not beta).
        int inc = (value == bestValue && ss->ply + 2 >= rootDepth && (int(nodes) & 14) == 0
                   && !is_win(std::abs(value) + 1));

        if (value + inc > bestValue)
        {
            bestValue = value;

            if (value + inc > alpha)
            {
                bestMove = move;

                if (PvNode && !rootNode)  // Update pv even in fail-high case
                    ss->pv->update(move, (ss + 1)->pv);

                if (value >= beta)
                {
                    // (*Scaler) Infrequent and small updates scale well
                    ss->cutoffCnt += (extension < 2) || PvNode;
                    assert(value >= beta);  // Fail high
                    break;
                }

                // Reduce other moves if we have found at least one score improvement
                if (depth > 2 && depth < 13 && !is_decisive(value))
                    depth -= 2;

                assert(depth > 0);
                alpha = value;  // Update alpha! Always alpha < beta
            }
        }

        // If the move is worse than some previously searched move,
        // remember it, to update its stats later.
        if (move != bestMove && moveCount <= SEARCHEDLIST_CAPACITY)
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

    if(CC && moveCount)
    {
	    bool T = bestValue > origAlpha;
	    bool T2 = V > origAlpha;
	    dbg_hit_on(T, 300);
	    dbg_hit_on(T2, 301);
	    dbg_hit_on(T, 310 + T2);
	    dbg_hit_on(T2, 320 + T);

	    dbg_hit_on(T, 400+C);
	    dbg_hit_on(T2, 410+C);
	    dbg_hit_on(T, 500 + 10*T2+C);
	    dbg_hit_on(T2, 520 + 10*T+C);

	    std::vector<bool>CL = {cutNode, improving, priorCapture, ss->inCheck, ttCapture, opponentWorsening, ss->ttHit, ss->ttPv,
	    (ss-1)->moveCount == 0, (ss-1)->currentMove == Move::null(), bool(excludedMove)};
	    dbg_correl_of(T, T2, 0);
	    dbg_correl_of(T, T2, origDepth);
	    
	    for(int i = 0; i < int(CL.size()); i++)
	    {
	         dbg_correl_of(T, T2, 1000*(i+1)+CL[i]*100);
	         dbg_correl_of(T, T2, 1000*(i+1)+CL[i]*100 + origDepth);
	    }
	    /*
	     * Hit #0: Total 2172719925 Hits 568475493 Hit Rate (%) 26.1642
Hit #1: Total 2172719925 Hits 579333760 Hit Rate (%) 26.664
Hit #2: Total 2172719925 Hits 179246322 Hit Rate (%) 8.24986
Hit #3: Total 2172719925 Hits 167447388 Hit Rate (%) 7.70681
Hit #4: Total 2172719925 Hits 139371885 Hit Rate (%) 6.41463
Hit #5: Total 2172719925 Hits 109744959 Hit Rate (%) 5.05104
Hit #6: Total 2172719925 Hits 85337951 Hit Rate (%) 3.9277
Hit #7: Total 2172719925 Hits 71422778 Hit Rate (%) 3.28725
Hit #8: Total 2172719925 Hits 59282293 Hit Rate (%) 2.72848
Hit #9: Total 2172719925 Hits 48310433 Hit Rate (%) 2.2235
Hit #10: Total 2172719925 Hits 38592026 Hit Rate (%) 1.77621
Hit #11: Total 2172719925 Hits 30520864 Hit Rate (%) 1.40473
Hit #12: Total 2172719925 Hits 23917146 Hit Rate (%) 1.10079
Hit #13: Total 2172719925 Hits 18468693 Hit Rate (%) 0.850026
Hit #14: Total 2172719925 Hits 14029544 Hit Rate (%) 0.645713
Hit #15: Total 2172719925 Hits 10494138 Hit Rate (%) 0.482995
Hit #16: Total 2172719925 Hits 7717682 Hit Rate (%) 0.355208
Hit #17: Total 2172719925 Hits 5632103 Hit Rate (%) 0.259219
Hit #18: Total 2172719925 Hits 4093430 Hit Rate (%) 0.188401
Hit #19: Total 2172719925 Hits 2983151 Hit Rate (%) 0.1373
Hit #20: Total 2172719925 Hits 2185987 Hit Rate (%) 0.100611
Hit #21: Total 2172719925 Hits 1614133 Hit Rate (%) 0.0742909
Hit #22: Total 2172719925 Hits 1202938 Hit Rate (%) 0.0553655
Hit #23: Total 2172719925 Hits 900776 Hit Rate (%) 0.0414584
Hit #24: Total 2172719925 Hits 680220 Hit Rate (%) 0.0313073
Hit #25: Total 2172719925 Hits 515997 Hit Rate (%) 0.0237489
Hit #26: Total 2172719925 Hits 386081 Hit Rate (%) 0.0177695
Hit #27: Total 2172719925 Hits 286443 Hit Rate (%) 0.0131836
Hit #28: Total 2172719925 Hits 207105 Hit Rate (%) 0.00953206
Hit #29: Total 2172719925 Hits 147442 Hit Rate (%) 0.00678606
Hit #30: Total 2172719925 Hits 101892 Hit Rate (%) 0.00468961
Hit #31: Total 2172719925 Hits 68872 Hit Rate (%) 0.00316985
Hit #101: Total 593917660 Hits 386247009 Hit Rate (%) 65.0338
Hit #102: Total 182488043 Hits 102913559 Hit Rate (%) 56.3947
Hit #103: Total 169289316 Hits 92765584 Hit Rate (%) 54.7971
Hit #104: Total 140530978 Hits 77123254 Hit Rate (%) 54.8799
Hit #105: Total 110451618 Hits 63019999 Hit Rate (%) 57.0567
Hit #106: Total 85832493 Hits 49100729 Hit Rate (%) 57.2053
Hit #107: Total 71830733 Hits 40573777 Hit Rate (%) 56.4853
Hit #108: Total 59631206 Hits 32467466 Hit Rate (%) 54.4471
Hit #109: Total 48581288 Hits 26220088 Hit Rate (%) 53.9716
Hit #110: Total 38617153 Hits 20632424 Hit Rate (%) 53.4281
Hit #111: Total 30521160 Hits 16173948 Hit Rate (%) 52.9926
Hit #112: Total 23917146 Hits 12567548 Hit Rate (%) 52.5462
Hit #113: Total 18468694 Hits 9630628 Hit Rate (%) 52.1457
Hit #114: Total 14029544 Hits 7219858 Hit Rate (%) 51.4618
Hit #115: Total 10494138 Hits 5278810 Hit Rate (%) 50.3025
Hit #116: Total 7717682 Hits 3764765 Hit Rate (%) 48.781
Hit #117: Total 5632103 Hits 2627196 Hit Rate (%) 46.6468
Hit #118: Total 4093430 Hits 1839166 Hit Rate (%) 44.9297
Hit #119: Total 2983151 Hits 1295343 Hit Rate (%) 43.422
Hit #120: Total 2185987 Hits 914976 Hit Rate (%) 41.8564
Hit #121: Total 1614133 Hits 655885 Hit Rate (%) 40.6339
Hit #122: Total 1202938 Hits 475931 Hit Rate (%) 39.5641
Hit #123: Total 900776 Hits 349825 Hit Rate (%) 38.836
Hit #124: Total 680220 Hits 255677 Hit Rate (%) 37.5874
Hit #125: Total 515997 Hits 188425 Hit Rate (%) 36.5167
Hit #126: Total 386081 Hits 133801 Hit Rate (%) 34.6562
Hit #127: Total 286443 Hits 94075 Hit Rate (%) 32.8425
Hit #128: Total 207105 Hits 62973 Hit Rate (%) 30.4063
Hit #129: Total 147442 Hits 41724 Hit Rate (%) 28.2986
Hit #130: Total 101892 Hits 25970 Hit Rate (%) 25.4878
Hit #131: Total 68872 Hits 16115 Hit Rate (%) 23.3985
Hit #201: Total 386247009 Hits 354373726 Hit Rate (%) 91.748
Hit #202: Total 102913559 Hits 91381489 Hit Rate (%) 88.7944
Hit #203: Total 92765584 Hits 83582696 Hit Rate (%) 90.101
Hit #204: Total 77123254 Hits 70210847 Hit Rate (%) 91.0372
Hit #205: Total 63019999 Hits 57881409 Hit Rate (%) 91.8461
Hit #206: Total 49100729 Hits 45469261 Hit Rate (%) 92.604
Hit #207: Total 40573777 Hits 37513059 Hit Rate (%) 92.4564
Hit #208: Total 32467466 Hits 29893785 Hit Rate (%) 92.073
Hit #209: Total 26220088 Hits 24090305 Hit Rate (%) 91.8773
Hit #210: Total 20632424 Hits 19128748 Hit Rate (%) 92.7121
Hit #211: Total 16173948 Hits 14989843 Hit Rate (%) 92.6789
Hit #212: Total 12567548 Hits 11611030 Hit Rate (%) 92.389
Hit #213: Total 9630628 Hits 8882294 Hit Rate (%) 92.2296
Hit #214: Total 7219858 Hits 6644602 Hit Rate (%) 92.0323
Hit #215: Total 5278810 Hits 4842992 Hit Rate (%) 91.744
Hit #216: Total 3764765 Hits 3441617 Hit Rate (%) 91.4165
Hit #217: Total 2627196 Hits 2386277 Hit Rate (%) 90.8298
Hit #218: Total 1839166 Hits 1656466 Hit Rate (%) 90.0661
Hit #219: Total 1295343 Hits 1154835 Hit Rate (%) 89.1528
Hit #220: Total 914976 Hits 803697 Hit Rate (%) 87.838
Hit #221: Total 655885 Hits 567752 Hit Rate (%) 86.5627
Hit #222: Total 475931 Hits 404557 Hit Rate (%) 85.0033
Hit #223: Total 349825 Hits 293590 Hit Rate (%) 83.9248
Hit #224: Total 255677 Hits 211634 Hit Rate (%) 82.774
Hit #225: Total 188425 Hits 154278 Hit Rate (%) 81.8777
Hit #226: Total 133801 Hits 108047 Hit Rate (%) 80.752
Hit #227: Total 94075 Hits 75408 Hit Rate (%) 80.1573
Hit #228: Total 62973 Hits 50054 Hit Rate (%) 79.4849
Hit #229: Total 41724 Hits 32839 Hit Rate (%) 78.7053
Hit #230: Total 25970 Hits 20432 Hit Rate (%) 78.6754
Hit #231: Total 16115 Hits 12614 Hit Rate (%) 78.2749
Hit #300: Total 23075010 Hits 4666807 Hit Rate (%) 20.2245
Hit #301: Total 23075010 Hits 2754497 Hit Rate (%) 11.9371
Hit #310: Total 20320513 Hits 2488410 Hit Rate (%) 12.2458
Hit #311: Total 2754497 Hits 2178397 Hit Rate (%) 79.0851
Hit #320: Total 18408203 Hits 576100 Hit Rate (%) 3.12958
Hit #321: Total 4666807 Hits 2178397 Hit Rate (%) 46.6785
Hit #400: Total 21868982 Hits 4017552 Hit Rate (%) 18.371
Hit #401: Total 1206028 Hits 649255 Hit Rate (%) 53.8342
Hit #410: Total 21868982 Hits 1548690 Hit Rate (%) 7.08167
Hit #411: Total 1206028 Hits 1205807 Hit Rate (%) 99.9817
Hit #500: Total 20320292 Hits 2488351 Hit Rate (%) 12.2456
Hit #501: Total 221 Hits 59 Hit Rate (%) 26.6968
Hit #510: Total 1548690 Hits 1529201 Hit Rate (%) 98.7416
Hit #511: Total 1205807 Hits 649196 Hit Rate (%) 53.8391
Hit #520: Total 17851430 Hits 19489 Hit Rate (%) 0.109173
Hit #521: Total 556773 Hits 556611 Hit Rate (%) 99.9709
Hit #530: Total 4017552 Hits 1529201 Hit Rate (%) 38.063
Hit #531: Total 649255 Hits 649196 Hit Rate (%) 99.9909
Correl. #0: Total 23075010 Coefficient 0.539517
Correl. #1: Total 14580182 Coefficient 0.689859
Correl. #2: Total 3241022 Coefficient 0.365404
Correl. #3: Total 1840848 Coefficient 0.311718
Correl. #4: Total 1158891 Coefficient 0.307967
Correl. #5: Total 706593 Coefficient 0.294102
Correl. #6: Total 494461 Coefficient 0.314058
Correl. #7: Total 407880 Coefficient 0.32036
Correl. #8: Total 348876 Coefficient 0.303959
Correl. #9: Total 270840 Coefficient 0.222567
Correl. #10: Total 25121 Coefficient 0.333996
Correl. #11: Total 295 Coefficient -nan
Correl. #13: Total 1 Coefficient -nan
Correl. #1000: Total 16021015 Coefficient 0.46684
Correl. #1001: Total 9937575 Coefficient 0.699869
Correl. #1002: Total 2475512 Coefficient 0.331821
Correl. #1003: Total 1114602 Coefficient 0.238875
Correl. #1004: Total 741799 Coefficient 0.180415
Correl. #1005: Total 490676 Coefficient 0.118077
Correl. #1006: Total 371212 Coefficient 0.114327
Correl. #1007: Total 325232 Coefficient 0.111671
Correl. #1008: Total 295059 Coefficient 0.108671
Correl. #1009: Total 246028 Coefficient 0.0950526
Correl. #1010: Total 23034 Coefficient 0.306817
Correl. #1011: Total 285 Coefficient -nan
Correl. #1013: Total 1 Coefficient -nan
Correl. #1100: Total 7053995 Coefficient 0.440862
Correl. #1101: Total 4642607 Coefficient 0.622772
Correl. #1102: Total 765510 Coefficient 0.229437
Correl. #1103: Total 726246 Coefficient 0.0687196
Correl. #1104: Total 417092 Coefficient 0.0304128
Correl. #1105: Total 215917 Coefficient -0.0506194
Correl. #1106: Total 123249 Coefficient -0.0293294
Correl. #1107: Total 82648 Coefficient 0.000636343
Correl. #1108: Total 53817 Coefficient 0.0243326
Correl. #1109: Total 24812 Coefficient 0.0574729
Correl. #1110: Total 2087 Coefficient 0.183125
Correl. #1111: Total 10 Coefficient -nan
Correl. #2000: Total 17586013 Coefficient 0.64321
Correl. #2001: Total 12282683 Coefficient 0.788415
Correl. #2002: Total 2364775 Coefficient 0.425842
Correl. #2003: Total 1056319 Coefficient 0.313085
Correl. #2004: Total 622271 Coefficient 0.282288
Correl. #2005: Total 387503 Coefficient 0.259637
Correl. #2006: Total 272740 Coefficient 0.271863
Correl. #2007: Total 226787 Coefficient 0.27834
Correl. #2008: Total 195156 Coefficient 0.263442
Correl. #2009: Total 159151 Coefficient 0.196706
Correl. #2010: Total 18401 Coefficient 0.340965
Correl. #2011: Total 227 Coefficient -nan
Correl. #2100: Total 5488997 Coefficient 0.317786
Correl. #2101: Total 2297499 Coefficient 0.373557
Correl. #2102: Total 876247 Coefficient 0.266116
Correl. #2103: Total 784529 Coefficient 0.283935
Correl. #2104: Total 536620 Coefficient 0.315803
Correl. #2105: Total 319090 Coefficient 0.306178
Correl. #2106: Total 221721 Coefficient 0.335119
Correl. #2107: Total 181093 Coefficient 0.343522
Correl. #2108: Total 153720 Coefficient 0.322586
Correl. #2109: Total 111689 Coefficient 0.235587
Correl. #2110: Total 6720 Coefficient 0.324687
Correl. #2111: Total 68 Coefficient -nan
Correl. #2113: Total 1 Coefficient -nan
Correl. #3000: Total 12240247 Coefficient 0.347622
Correl. #3001: Total 5536104 Coefficient 0.455906
Correl. #3002: Total 1871832 Coefficient 0.27064
Correl. #3003: Total 1608304 Coefficient 0.277223
Correl. #3004: Total 1094417 Coefficient 0.295828
Correl. #3005: Total 672246 Coefficient 0.284986
Correl. #3006: Total 468275 Coefficient 0.306759
Correl. #3007: Total 386893 Coefficient 0.31598
Correl. #3008: Total 330378 Coefficient 0.301141
Correl. #3009: Total 255421 Coefficient 0.221451
Correl. #3010: Total 16206 Coefficient 0.326005
Correl. #3011: Total 170 Coefficient -nan
Correl. #3013: Total 1 Coefficient -nan
Correl. #3100: Total 10834763 Coefficient 0.922228
Correl. #3101: Total 9044078 Coefficient 0.944374
Correl. #3102: Total 1369190 Coefficient 0.722878
Correl. #3103: Total 232544 Coefficient 0.421015
Correl. #3104: Total 64474 Coefficient 0.454877
Correl. #3105: Total 34347 Coefficient 0.462671
Correl. #3106: Total 26186 Coefficient 0.451762
Correl. #3107: Total 20987 Coefficient 0.37106
Correl. #3108: Total 18498 Coefficient 0.33574
Correl. #3109: Total 15419 Coefficient 0.212763
Correl. #3110: Total 8915 Coefficient 0.391019
Correl. #3111: Total 125 Coefficient -nan
Correl. #4000: Total 23075010 Coefficient 0.539517
Correl. #4001: Total 14580182 Coefficient 0.689859
Correl. #4002: Total 3241022 Coefficient 0.365404
Correl. #4003: Total 1840848 Coefficient 0.311718
Correl. #4004: Total 1158891 Coefficient 0.307967
Correl. #4005: Total 706593 Coefficient 0.294102
Correl. #4006: Total 494461 Coefficient 0.314058
Correl. #4007: Total 407880 Coefficient 0.32036
Correl. #4008: Total 348876 Coefficient 0.303959
Correl. #4009: Total 270840 Coefficient 0.222567
Correl. #4010: Total 25121 Coefficient 0.333996
Correl. #4011: Total 295 Coefficient -nan
Correl. #4013: Total 1 Coefficient -nan
Correl. #5000: Total 22280927 Coefficient 0.552445
Correl. #5001: Total 14297276 Coefficient 0.691306
Correl. #5002: Total 3092839 Coefficient 0.363363
Correl. #5003: Total 1687303 Coefficient 0.331122
Correl. #5004: Total 1072725 Coefficient 0.329869
Correl. #5005: Total 658511 Coefficient 0.309737
Correl. #5006: Total 465651 Coefficient 0.322424
Correl. #5007: Total 386333 Coefficient 0.322235
Correl. #5008: Total 333361 Coefficient 0.29642
Correl. #5009: Total 261929 Coefficient 0.20588
Correl. #5010: Total 24703 Coefficient 0.322414
Correl. #5011: Total 295 Coefficient -nan
Correl. #5013: Total 1 Coefficient -nan
Correl. #5100: Total 794083 Coefficient 0.179993
Correl. #5101: Total 282906 Coefficient 0.551237
Correl. #5102: Total 148183 Coefficient 0.113877
Correl. #5103: Total 153545 Coefficient -0.0256848
Correl. #5104: Total 86166 Coefficient -0.0126726
Correl. #5105: Total 48082 Coefficient 0.0405643
Correl. #5106: Total 28810 Coefficient 0.117872
Correl. #5107: Total 21547 Coefficient 0.187375
Correl. #5108: Total 15515 Coefficient 0.257964
Correl. #5109: Total 8911 Coefficient 0.295004
Correl. #5110: Total 418 Coefficient 0.239049
Correl. #6000: Total 18068077 Coefficient 0.653158
Correl. #6001: Total 12636876 Coefficient 0.776082
Correl. #6002: Total 2400031 Coefficient 0.400757
Correl. #6003: Total 1059185 Coefficient 0.369601
Correl. #6004: Total 613100 Coefficient 0.402403
Correl. #6005: Total 375282 Coefficient 0.409585
Correl. #6006: Total 292159 Coefficient 0.383385
Correl. #6007: Total 256040 Coefficient 0.36413
Correl. #6008: Total 231127 Coefficient 0.303306
Correl. #6009: Total 190750 Coefficient 0.205758
Correl. #6010: Total 13262 Coefficient 0.307755
Correl. #6011: Total 265 Coefficient -nan
Correl. #6100: Total 5006933 Coefficient 0.267736
Correl. #6101: Total 1943306 Coefficient 0.389676
Correl. #6102: Total 840991 Coefficient 0.283984
Correl. #6103: Total 781663 Coefficient 0.231842
Correl. #6104: Total 545791 Coefficient 0.207467
Correl. #6105: Total 331311 Coefficient 0.172537
Correl. #6106: Total 202302 Coefficient 0.217714
Correl. #6107: Total 151840 Coefficient 0.247042
Correl. #6108: Total 117749 Coefficient 0.275283
Correl. #6109: Total 80090 Coefficient 0.221541
Correl. #6110: Total 11859 Coefficient 0.337655
Correl. #6111: Total 30 Coefficient -nan
Correl. #6113: Total 1 Coefficient -nan
Correl. #7100: Total 23075010 Coefficient 0.539517
Correl. #7101: Total 14580182 Coefficient 0.689859
Correl. #7102: Total 3241022 Coefficient 0.365404
Correl. #7103: Total 1840848 Coefficient 0.311718
Correl. #7104: Total 1158891 Coefficient 0.307967
Correl. #7105: Total 706593 Coefficient 0.294102
Correl. #7106: Total 494461 Coefficient 0.314058
Correl. #7107: Total 407880 Coefficient 0.32036
Correl. #7108: Total 348876 Coefficient 0.303959
Correl. #7109: Total 270840 Coefficient 0.222567
Correl. #7110: Total 25121 Coefficient 0.333996
Correl. #7111: Total 295 Coefficient -nan
Correl. #7113: Total 1 Coefficient -nan
Correl. #8000: Total 20519073 Coefficient 0.599843
Correl. #8001: Total 13975411 Coefficient 0.712151
Correl. #8002: Total 2872798 Coefficient 0.413404
Correl. #8003: Total 1306429 Coefficient 0.365306
Correl. #8004: Total 800352 Coefficient 0.340497
Correl. #8005: Total 491003 Coefficient 0.322426
Correl. #8006: Total 345672 Coefficient 0.337759
Correl. #8007: Total 283411 Coefficient 0.339441
Correl. #8008: Total 237903 Coefficient 0.295599
Correl. #8009: Total 182375 Coefficient 0.214104
Correl. #8010: Total 23424 Coefficient 0.361188
Correl. #8011: Total 294 Coefficient -nan
Correl. #8013: Total 1 Coefficient -nan
Correl. #8100: Total 2555937 Coefficient 0.227938
Correl. #8101: Total 604771 Coefficient 0.325766
Correl. #8102: Total 368224 Coefficient 0.193533
Correl. #8103: Total 534419 Coefficient 0.125398
Correl. #8104: Total 358539 Coefficient 0.1907
Correl. #8105: Total 215590 Coefficient 0.235811
Correl. #8106: Total 148789 Coefficient 0.262657
Correl. #8107: Total 124469 Coefficient 0.287611
Correl. #8108: Total 110973 Coefficient 0.299476
Correl. #8109: Total 88465 Coefficient 0.231961
Correl. #8110: Total 1697 Coefficient 0.189095
Correl. #8111: Total 1 Coefficient -nan
Correl. #9000: Total 20639998 Coefficient 0.532508
Correl. #9001: Total 12560106 Coefficient 0.686727
Correl. #9002: Total 2869252 Coefficient 0.357917
Correl. #9003: Total 1808953 Coefficient 0.309413
Correl. #9004: Total 1155116 Coefficient 0.30801
Correl. #9005: Total 702465 Coefficient 0.293991
Correl. #9006: Total 493595 Coefficient 0.314062
Correl. #9007: Total 407411 Coefficient 0.320497
Correl. #9008: Total 347728 Coefficient 0.304414
Correl. #9009: Total 270052 Coefficient 0.222881
Correl. #9010: Total 25024 Coefficient 0.334182
Correl. #9011: Total 295 Coefficient -nan
Correl. #9013: Total 1 Coefficient -nan
Correl. #9100: Total 2435012 Coefficient 0.18365
Correl. #9101: Total 2020076 Coefficient 0.233551
Correl. #9102: Total 371770 Coefficient 0.125026
Correl. #9103: Total 31895 Coefficient 0.143392
Correl. #9104: Total 3775 Coefficient 0.146611
Correl. #9105: Total 4128 Coefficient 0.184651
Correl. #9106: Total 866 Coefficient 0.19392
Correl. #9107: Total 469 Coefficient 0.120708
Correl. #9108: Total 1148 Coefficient -0.0201654
Correl. #9109: Total 788 Coefficient -nan
Correl. #9110: Total 97 Coefficient -nan
Correl. #10000: Total 23003724 Coefficient 0.539868
Correl. #10001: Total 14526957 Coefficient 0.690257
Correl. #10002: Total 3233982 Coefficient 0.366088
Correl. #10003: Total 1834982 Coefficient 0.312366
Correl. #10004: Total 1157901 Coefficient 0.308114
Correl. #10005: Total 703019 Coefficient 0.293926
Correl. #10006: Total 494012 Coefficient 0.313986
Correl. #10007: Total 407780 Coefficient 0.320358
Correl. #10008: Total 348856 Coefficient 0.303958
Correl. #10009: Total 270818 Coefficient 0.222564
Correl. #10010: Total 25121 Coefficient 0.333996
Correl. #10011: Total 295 Coefficient -nan
Correl. #10013: Total 1 Coefficient -nan
Correl. #10100: Total 71286 Coefficient 0.326172
Correl. #10101: Total 53225 Coefficient 0.452385
Correl. #10102: Total 7040 Coefficient 0.0666694
Correl. #10103: Total 5866 Coefficient 0.0902175
Correl. #10104: Total 990 Coefficient 0.0344282
Correl. #10105: Total 3574 Coefficient 0.218642
Correl. #10106: Total 449 Coefficient 0.160024
Correl. #10107: Total 100 Coefficient -nan
Correl. #10108: Total 20 Coefficient -nan
Correl. #10109: Total 22 Coefficient -nan
Correl. #11000: Total 21868982 Coefficient 0.572962
Correl. #11001: Total 14580182 Coefficient 0.689859
Correl. #11002: Total 3083440 Coefficient 0.380795
Correl. #11003: Total 1381474 Coefficient 0.210163
Correl. #11004: Total 868735 Coefficient 0.187726
Correl. #11005: Total 549902 Coefficient 0.174152
Correl. #11006: Total 416867 Coefficient 0.165342
Correl. #11007: Total 365994 Coefficient 0.157451
Correl. #11008: Total 329539 Coefficient 0.151391
Correl. #11009: Total 267456 Coefficient 0.151947
Correl. #11010: Total 25097 Coefficient 0.347355
Correl. #11011: Total 295 Coefficient -nan
Correl. #11013: Total 1 Coefficient -nan
Correl. #11100: Total 1206028 Coefficient 0.00736945
Correl. #11102: Total 157582 Coefficient 0.0143204
Correl. #11103: Total 459374 Coefficient 0.00418267
Correl. #11104: Total 290156 Coefficient -nan
Correl. #11105: Total 156691 Coefficient -nan
Correl. #11106: Total 77594 Coefficient -nan
Correl. #11107: Total 41886 Coefficient -nan
Correl. #11108: Total 19337 Coefficient -nan
Correl. #11109: Total 3384 Coefficient -nan
Correl. #11110: Total 24 Coefficient -nan

===========================
Total time (ms) : 4118828
Nodes searched  : 1991860457
Nodes/second    : 483598
	     */
    }

    // Adjust best value for fail high cases
    if (bestValue >= beta && !is_decisive(bestValue) && !is_decisive(alpha))
        bestValue = (bestValue * depth + beta) / (depth + 1);

    if (!moveCount)
        bestValue = excludedMove ? alpha : ss->inCheck ? mated_in(ss->ply) : VALUE_DRAW;

    // If there is a move that produces search value greater than alpha,
    // we update the stats of searched moves.
    else if (bestMove)
    {
        update_all_stats(pos, ss, *this, bestMove, prevSq, quietsSearched, capturesSearched, depth,
                         ttData.move, PvNode);
        if (!PvNode)
            ttMoveHistory << (bestMove == ttData.move ? 792 : -779);
    }

    // Bonus for prior quiet countermove that caused the fail low
    else if (!priorCapture && prevSq != SQ_NONE)
    {
        int bonusScale = -245;
        bonusScale -= (ss - 1)->statScore / 98;
        bonusScale += std::min(59 * depth, 430);
        bonusScale += 191 * ((ss - 1)->moveCount > 8);
        bonusScale += 143 * (!ss->inCheck && bestValue <= ss->staticEval - 103);
        bonusScale += 151 * (!(ss - 1)->inCheck && bestValue <= -(ss - 1)->staticEval - 78);

        bonusScale = std::max(bonusScale, 0);

        // scaledBonus ranges from 0 to roughly 2.3M, overflows happen for multipliers larger than 900
        const int scaledBonus = std::min(141 * depth - 82, 1472) * bonusScale;

        update_continuation_histories(ss - 1, pos.piece_on(prevSq), prevSq,
                                      scaledBonus * 236 / 16384);

        mainHistory[~us][((ss - 1)->currentMove).raw()] << scaledBonus * 234 / 32768;

        if (type_of(pos.piece_on(prevSq)) != PAWN && ((ss - 1)->currentMove).type_of() != PROMOTION)
            sharedHistory.pawn_entry(pos)[pos.piece_on(prevSq)][prevSq] << scaledBonus * 322 / 8192;
    }

    // Bonus for prior capture countermove that caused the fail low
    else if (priorCapture && prevSq != SQ_NONE)
    {
        Piece capturedPiece = pos.captured_piece();
        assert(capturedPiece != NO_PIECE);
        captureHistory[pos.piece_on(prevSq)][prevSq][type_of(capturedPiece)] << 901;
    }

    if (PvNode)
        bestValue = std::min(bestValue, maxValue);

    // If no good move is found and the previous position was ttPv, then the previous
    // opponent move is probably good and the new position is added to the search tree.
    if (bestValue <= alpha)
        ss->ttPv = ss->ttPv || (ss - 1)->ttPv;

    // Write gathered information in transposition table. Note that the
    // static evaluation is saved as it was before correction history.
    if (!excludedMove && !(rootNode && pvIdx))
        ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), ss->ttPv,
                       bestValue >= beta    ? BOUND_LOWER
                       : PvNode && bestMove ? BOUND_EXACT
                                            : BOUND_UPPER,
                       moveCount != 0 ? depth : std::min(MAX_PLY - 1, depth + 6), bestMove,
                       unadjustedStaticEval, tt.generation());

    // Adjust correction history if the best move is not a capture
    // and the error direction matches whether we are above/below bounds.
    if (!ss->inCheck && !(bestMove && pos.capture(bestMove))
        && (bestValue > ss->staticEval) == bool(bestMove))
    {
        auto bonus =
          std::clamp(int(bestValue - ss->staticEval) * depth * (bestMove ? 12 : 18) / 128,
                     -CORRECTION_HISTORY_LIMIT / 4, CORRECTION_HISTORY_LIMIT / 4);
        update_correction_history(pos, ss, *this, 1114 * bonus / 1024);
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

    for(int d = 0; d < 32; d++)
         dbg_hit_on(d == 0, d);

    assert(alpha >= -VALUE_INFINITE && alpha < beta && beta <= VALUE_INFINITE);
    assert(PvNode || (alpha == beta - 1));

    // Check if we have an upcoming move that draws by repetition
    if (alpha < VALUE_DRAW && pos.upcoming_repetition(ss->ply))
    {
        alpha = value_draw(nodes);
        if (alpha >= beta)
            return alpha;
    }

    PVMoves   pv;
    StateInfo st;

    Key   posKey;
    Move  move, bestMove;
    Value bestValue, value, futilityBase;
    bool  pvHit, givesCheck, capture;
    int   moveCount;

    // Step 1. Initialize node
    if (PvNode)
    {
        (ss + 1)->pv = &pv;
        ss->pv->clear();
    }

    bestMove    = Move::none();
    ss->inCheck = pos.checkers();
    moveCount   = 0;

    // Used to send selDepth info to GUI (selDepth counts from 1, ply from 0)
    if (PvNode && selDepth < ss->ply + 1)
        selDepth = ss->ply + 1;

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
    Value unadjustedStaticEval = VALUE_NONE;
    if (ss->inCheck)
        bestValue = futilityBase = -VALUE_INFINITE;
    else
    {
        const auto correctionValue = correction_value(*this, pos, ss);

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
            unadjustedStaticEval = evaluate(pos);
            ss->staticEval       = bestValue =
              to_corrected_static_eval(unadjustedStaticEval, correctionValue);
        }

        // Stand pat. Return immediately if static value is at least beta
        if (bestValue >= beta)
        {
            if (!is_decisive(bestValue))
                bestValue = (467 * bestValue + 557 * beta) / 1024;

            if (!ss->ttHit)
                ttWriter.write(posKey, VALUE_NONE, false, BOUND_LOWER, DEPTH_UNSEARCHED,
                               Move::none(), unadjustedStaticEval, tt.generation());
            return bestValue;
        }

        if (bestValue > alpha)
            alpha = bestValue;

        futilityBase = ss->staticEval + 335;
    }

    const PieceToHistory* contHist[] = {(ss - 1)->continuationHistory};

    Square prevSq = ((ss - 1)->currentMove).is_ok() ? ((ss - 1)->currentMove).to_sq() : SQ_NONE;

    // Initialize a MovePicker object for the current position, and prepare to search
    // the moves. We presently use two stages of move generator in quiescence search:
    // captures, or evasions only when in check.
    MovePicker mp(pos, ttData.move, DEPTH_QS, &mainHistory, &lowPlyHistory, &captureHistory,
                  contHist, &sharedHistory, ss->ply);

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
                    bestValue = std::max(bestValue, std::min(alpha, futilityBase));
                    continue;
                }
            }

            // Skip non-captures
            if (!capture)
                continue;

            // Do not search moves with bad enough SEE values
            if (!pos.see_ge(move, -74))
                continue;
        }

        // Step 7. Make and search the move
        do_move(pos, move, st, givesCheck, ss);

        value = -qsearch<nodeType>(pos, ss + 1, -beta, -alpha);
        undo_move(pos, move);

        assert(value > -VALUE_INFINITE && value < VALUE_INFINITE);

        // Step 8. Check for a new best move
        if (value > bestValue)
        {
            bestValue = value;

            if (value > alpha)
            {
                bestMove = move;

                if (PvNode)  // Update pv even in fail-high case
                    ss->pv->update(move, (ss + 1)->pv);

                if (value < beta)  // Update alpha here!
                    alpha = value;
                else
                    break;  // Fail high
            }
        }
    }

    // Step 9. Check for mate and stalemate
    // All legal moves have been searched. A special case: if we are
    // in check and no legal moves were found, it is checkmate.
    if (!moveCount)
    {
        if (ss->inCheck)  // Checkmate!
        {
            assert(!MoveList<LEGAL>(pos).size());
            return mated_in(ss->ply);  // Plies to mate from the root
        }

        // Only check for stalemate under specific conditions
        Color us = pos.side_to_move();
        if (!(pawn_single_push_bb(us, pos.pieces(us, PAWN)) & ~pos.pieces())
            && !pos.non_pawn_material(us) && type_of(pos.captured_piece()) >= KNIGHT
            && !MoveList<LEGAL>(pos).size())
            bestValue = VALUE_DRAW;
    }

    if (!is_decisive(bestValue) && bestValue > beta)
        bestValue = (481 * bestValue + 543 * beta) / 1024;

    // Save gathered info in transposition table. The static evaluation
    // is saved as it was before adjustment by correction history.
    ttWriter.write(posKey, value_to_tt(bestValue, ss->ply), pvHit,
                   bestValue >= beta ? BOUND_LOWER : BOUND_UPPER, DEPTH_QS, bestMove,
                   unadjustedStaticEval, tt.generation());

    assert(bestValue > -VALUE_INFINITE && bestValue < VALUE_INFINITE);

    return bestValue;
}

int Search::Worker::reduction(bool i, Depth d, int mn, int delta) const {
    int reductionScale = reductions[d] * reductions[mn];
    return reductionScale - delta * 617 / rootDelta + !i * reductionScale * 194 / 512 + 1027;
}

// elapsed() returns the time elapsed since the search started. If the
// 'nodestime' option is enabled, it will return the count of nodes searched
// instead. This function is called to check whether the search should be
// stopped based on predefined thresholds like time limits or nodes searched.
TimePoint Search::Worker::elapsed() const {
    return main_manager()->tm.elapsed([this]() { return threads.nodes_searched(); });
}

Value Search::Worker::evaluate(const Position& pos) {
    return Eval::evaluate(network[numaAccessToken], pos, accumulatorStack, refreshTable,
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
        if (is_mate(v) && VALUE_MATE - v > 100 - r50c)
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
        if (is_mated(v) && VALUE_MATE + v > 100 - r50c)
            return VALUE_TB_LOSS_IN_MAX_PLY + 1;

        // Downgrade a potentially false TB score.
        if (VALUE_TB + v > 100 - r50c)
            return VALUE_TB_LOSS_IN_MAX_PLY + 1;

        return v + ply;
    }

    return v;
}


// Updates stats at the end of search() when a bestMove is found
void update_all_stats(const Position& pos,
                      Stack*          ss,
                      Search::Worker& workerThread,
                      Move            bestMove,
                      Square          prevSq,
                      SearchedList&   quietsSearched,
                      SearchedList&   capturesSearched,
                      Depth           depth,
                      Move            ttMove,
                      bool            PvNode) {

    CapturePieceToHistory& captureHistory = workerThread.captureHistory;
    Piece                  movedPiece     = pos.moved_piece(bestMove);
    PieceType              capturedPiece;

    int bonus =
      std::min(134 * depth - 79, 1572) + 382 * (bestMove == ttMove) + (ss - 1)->statScore / 30;
    int malus = std::min(1005 * depth - 205, 2218);

    if (!PvNode)
        // Important: don't remove the cast to a 64-bit number else the multiplication
        // can overflow on 32-bit platforms which would change the bench signature
        bonus += int(bonus * u64(quietsSearched.size() + capturesSearched.size()) / 256);

    if (!pos.capture_stage(bestMove))
    {
        update_quiet_histories(pos, ss, workerThread, bestMove, bonus * 824 / 1024);

        int actualMalus = malus * 1136 / 1024;
        // Decrease stats for all non-best quiet moves
        for (Move move : quietsSearched)
        {
            actualMalus = actualMalus * 956 / 1024;
            update_quiet_histories(pos, ss, workerThread, move, -actualMalus);
        }
    }
    else
    {
        // Increase stats for the best move in case it was a capture move
        capturedPiece = type_of(pos.piece_on(bestMove.to_sq()));
        captureHistory[movedPiece][bestMove.to_sq()][capturedPiece] << bonus * 1366 / 1024;
    }

    // Extra penalty for a quiet early move that was not a TT move in
    // previous ply when it gets refuted.
    if (prevSq != SQ_NONE && ((ss - 1)->moveCount == 1 + (ss - 1)->ttHit) && !pos.captured_piece())
        update_continuation_histories(ss - 1, pos.piece_on(prevSq), prevSq, -malus * 683 / 1024);

    // Decrease stats for all non-best capture moves
    for (Move move : capturesSearched)
    {
        movedPiece    = pos.moved_piece(move);
        capturedPiece = type_of(pos.piece_on(move.to_sq()));
        captureHistory[movedPiece][move.to_sq()][capturedPiece] << -malus * 1518 / 1024;
    }
}


// Updates the continuation histories for the move pairs formed by
// the current move and the moves played in previous plies.
void update_continuation_histories(Stack* ss, Piece pc, Square to, int bonus) {
    static constexpr std::array<ConthistBonus, 6> conthist_bonuses = {
      {{1, 1040}, {2, 780}, {3, 300}, {4, 537}, {5, 129}, {6, 423}}};

    // Multipliers for positive history consistency
    constexpr int CMHCMultipliers[] = {96, 113, 101, 105, 127, 121, 126};
    int           positiveCount     = 0;

    for (const auto [i, weight] : conthist_bonuses)
    {
        // Only update the first 2 continuation histories if we are in check
        if (ss->inCheck && i > 2)
            break;

        if (((ss - i)->currentMove).is_ok())
        {
            auto& historyEntry = (*(ss - i)->continuationHistory)[pc][to];
            if (historyEntry > 0)
                positiveCount++;

            int multiplier = CMHCMultipliers[positiveCount];
            historyEntry << (bonus * weight * multiplier / 131072) + 71 * (i < 2);
        }
    }
}

// Updates move sorting heuristics

void update_quiet_histories(
  const Position& pos, Stack* ss, Search::Worker& workerThread, Move move, int bonus) {

    Color us = pos.side_to_move();
    workerThread.mainHistory[us][move.raw()] << bonus;  // Untuned to prevent duplicate effort

    if (ss->ply < LOW_PLY_HISTORY_SIZE)
        workerThread.lowPlyHistory[ss->ply][move.raw()] << bonus * 663 / 1024;

    update_continuation_histories(ss, pos.moved_piece(move), move.to_sq(), bonus * 820 / 1024);

    workerThread.sharedHistory.pawn_entry(pos)[pos.moved_piece(move)][move.to_sq()]
      << bonus * (bonus > -7 ? 1038 : 525) / 1024;
}
}

// When playing with strength handicap, choose the best move among a set of
// RootMoves using a statistical rule dependent on 'level'. Idea by Heinz van Saanen.
Move Skill::pick_best(const RootMoves& rootMoves, usize multiPV) {
    static PRNG rng(now());  // PRNG sequence should be non-deterministic

    // With tablebases at the root, rootMoves are ordered by tbRank rather than by
    // score, so compute the score range explicitly to keep 'delta' non-negative.
    Value topScore = rootMoves[0].score;
    Value minScore = rootMoves[0].score;
    for (usize i = 1; i < multiPV; ++i)
    {
        topScore = std::max(topScore, rootMoves[i].score);
        minScore = std::min(minScore, rootMoves[i].score);
    }
    int    delta    = std::min(topScore - minScore, int(PawnValue));
    int    maxScore = -VALUE_INFINITE;
    double weakness = 120 - 2 * level;

    // Choose best move. For each move score we add two terms, both dependent on
    // weakness. One is deterministic and bigger for weaker levels, and one is
    // random. Then we choose the move with the resulting highest score.
    for (usize i = 0; i < multiPV; ++i)
    {
        // This is our magic formula
        int push = int(weakness * int(topScore - rootMoves[i].score)
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

    if ((worker.limits.use_time_management() && (elapsed > tm.maximum() || stopOnPonderhit))
        || (worker.limits.movetime && elapsed >= worker.limits.movetime)
        || (worker.limits.nodes && worker.threads.nodes_searched() >= worker.limits.nodes))
        worker.threads.stop = true;
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
    while (usize(ply) < rootMove.pv.size())
    {
        Move& pvMove = rootMove.pv[ply];

        RootMoves legalMoves;
        for (const auto& m : MoveList<LEGAL>(pos))
            legalMoves.emplace_back(m);

        TB::Config config = TB::rank_root_moves(options, pos, legalMoves, false, time_abort);
        RootMove&  rm     = *std::find(legalMoves.begin(), legalMoves.end(), pvMove);

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
        TB::Config config = TB::rank_root_moves(options, pos, legalMoves, true, time_abort);

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
    // rule if it has been reached on the board with a non-optimal 50 move counter
    // (e.g. 8/8/6k1/3B4/3K4/4N3/8/8 w - - 54 106 ) which TB with dtz counter rounding
    // cannot always correctly rank. See also
    // https://github.com/official-stockfish/Stockfish/issues/5175#issuecomment-2058893495
    // We adjust the score to match the found PV. Note that a TB loss score can be
    // displayed if the engine did not find a drawing move yet, but eventually search
    // will figure it out (e.g. 1kq5/q2r4/5K2/8/8/8/8/7Q w - - 96 1 )
    if (pos.is_draw(0))
        v = VALUE_DRAW;

    // Undo the PV moves
    for (usize i = rootMove.pv.size(); i > 0; --i)
        pos.undo_move(rootMove.pv[i - 1]);

    // Inform if we couldn't get a full extension in time
    if (time_abort())
        sync_cout
          << "info string Syzygy based PV extension requires more time, increase Move Overhead as needed."
          << sync_endl;
}

void SearchManager::output_pv(Search::Worker&           worker,
                              const ThreadPool&         threads,
                              const TranspositionTable& tt,
                              Depth                     depth) {

    const auto nodes     = threads.nodes_searched();
    auto&      rootMoves = worker.rootMoves;
    auto&      pos       = worker.rootPos;
    usize      multiPV   = std::min(usize(worker.options["MultiPV"]), rootMoves.size());
    u64        tbHits    = threads.tb_hits() + (worker.tbConfig.rootInTB ? rootMoves.size() : 0);

    for (usize i = 0; i < multiPV; ++i)
    {
        bool usePreviousScore = rootMoves[i].score == -VALUE_INFINITE;

        if (depth == 1 && usePreviousScore && i > 0)
            continue;

        Depth d = usePreviousScore ? std::max(1, depth - 1) : depth;
        Value v = usePreviousScore ? rootMoves[i].previousScore : rootMoves[i].uciScore;

        if (v == -VALUE_INFINITE)
            v = VALUE_ZERO;

        bool isTBScore = worker.tbConfig.rootInTB && !is_mate_or_mated(v);
        v              = isTBScore ? rootMoves[i].tbScore : v;

        // Potentially correct and extend the PV, and in exceptional cases v.
        // Previous PVs have already been extended. Bound flags indicate an unreliable PV.
        if (is_decisive(v) && !is_mate_or_mated(v) && !usePreviousScore
            && (!rootMoves[i].score_is_bound() || isTBScore))
            syzygy_extend_pv(worker.options, worker.limits, pos, rootMoves[i], v);

        std::string pv;
        for (Move m : usePreviousScore ? rootMoves[i].previousPV : rootMoves[i].pv)
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

        // TB and previous scores are exact, even though their bound flags may say otherwise.
        if (!(isTBScore || usePreviousScore))
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

    assert(pv.size() == 1 && pv[0] != Move::none());

    StateInfo st;
    pos.do_move(pv[0], st, &tt);

    if (!pos.is_draw(1))
    {
        auto [ttHit, ttData, ttWriter] = tt.probe(pos.key());
        if (ttHit && MoveList<LEGAL>(pos).contains(ttData.move))
            pv.push_back(ttData.move);
    }

    pos.undo_move(pv[0]);
    return pv.size() > 1;
}


}  // namespace Stockfish
