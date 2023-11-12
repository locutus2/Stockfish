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

#ifndef STATS_H_INCLUDED
#define STATS_H_INCLUDED

#include <tuple>
#include <vector>

namespace Stockfish {

enum HistoryType : int {
    HISTORY_MAIN,
    HISTORY_PAWN,
    HISTORY_INCHECK,
    HISTORY_CMH0,
    HISTORY_CMH1,
    HISTORY_CMH2,
    HISTORY_CMH3,
    HISTORY_CMH0_POS,
    HISTORY_CMH0_NEG,
    HISTORY_MAIN_PAWN,
    HISTORY_MAIN_PAWN_SHIFT,
    HISTORY_MAIN_SHIFT_PAWN_SHIFT,
    N_HISTORY
};

extern std::vector<int> HISTORY_SCALE;
extern std::vector<int> HISTORY_WEIGHT;

extern int Dmax;
extern int Dmin;

constexpr int HISTORY_DIVISOR[N_HISTORY]      = {7183,  8192,  7183,  29952, 29952, 29952,
                                                 29952, 14976, 14976, 7183,  7183,  7183};
constexpr int HISTORY_SCALE_START[N_HISTORY]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
constexpr int HISTORY_WEIGHT_START[N_HISTORY] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

void init_stats(bool onlyD = false);


//---------------
// search + measurement

constexpr int HISTORY_BUCKETS = 1000;

constexpr bool STATS_EVASION_MAIN = true;
constexpr bool STATS_EVASION_QS   = false;
constexpr bool USE_DEPTH_WEIGHT   = true;


//---------------
// uci stats command

constexpr std::tuple<int, int, const char*> STATS_STEPS[] = {
  {-2, 1, "-2"},  {-1, 1, "-1"}, {-1, 2, "-0.5"}, {-1, 4, "-0.25"}, {0, 1, "0"},
  {1, 4, "0.25"}, {1, 2, "0.5"}, {1, 1, "1"},     {2, 1, "2"},
};

constexpr std::tuple<int, const char*> STATS_PARAMS[] = {
  {HISTORY_MAIN, "main"}, {HISTORY_PAWN, "pawn"}, {HISTORY_INCHECK, "incheck"},
  {HISTORY_CMH0, "cmh0"}, {HISTORY_CMH1, "cmh1"}, {HISTORY_CMH2, "cmh2"},
  {HISTORY_CMH3, "cmh3"},
  //{HISTORY_CMH0_POS, "cmh0_pos"},
  //{HISTORY_CMH0_NEG, "cmh0_neg"},
};

}
#endif
