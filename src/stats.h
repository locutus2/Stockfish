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

#include <algorithm>
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
    HISTORY_CMH4,
    HISTORY_CMH5,
    HISTORY_CAPTURE,
    HISTORY_CMH0_POS,
    HISTORY_CMH0_NEG,
    HISTORY_MAIN_PAWN,
    HISTORY_MAIN_PAWN_SHIFT,
    HISTORY_MAIN_SHIFT_PAWN_SHIFT,
    HISTORY_REF_ORDER,
    HISTORY_GIVES_CHECK,
    HISTORY_ESCAPE,
    HISTORY_EN_PRISE,
    N_HISTORY
};

extern std::vector<int> HISTORY_SCALE;
extern std::vector<int> HISTORY_WEIGHT;

extern int Dmax;
extern int Dmin;

constexpr int HISTORY_DIVISOR[N_HISTORY] = {7183,  8192,     7183,      29952, 29952, 29952, 29952,
                                            29952, 29952,    10692,     14976, 14976, 7183,  7183,
                                            7183,  7183 / 8, 16384 / 2, 25000, 25000};

constexpr int HISTORY_SCALE_QUIET_MASTER[N_HISTORY]  = {1, 1, 1, 1, 1, 4, 1, 1, 1, 1,
                                                        1, 1, 1, 1, 1, 1, 1, 1, 1};
constexpr int HISTORY_WEIGHT_QUIET_MASTER[N_HISTORY] = {2, 2, 0, 2, 1, 1, 1, 0, 1, 0,
                                                        0, 0, 0, 0, 0, 0, 1, 1, -1};

constexpr int HISTORY_SCALE_QUIET_EVASION_MASTER[N_HISTORY]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                                1, 1, 1, 1, 1, 1, 1, 1, 1};
constexpr int HISTORY_WEIGHT_QUIET_EVASION_MASTER[N_HISTORY] = {1, 1, 0, 1, 0, 0, 0, 0, 0, 0,
                                                                0, 0, 0, 0, 0, 0, 0, 0, 0};

constexpr int HISTORY_SCALE_CAPTURE_EVASION_MASTER[N_HISTORY]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                                  1, 1, 1, 1, 1, 1, 1, 1, 1};
constexpr int HISTORY_WEIGHT_CAPTURE_EVASION_MASTER[N_HISTORY] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                                  0, 0, 0, 0, 0, 0, 0, 0, 0};

constexpr int HISTORY_SCALE_REFUTATION_MASTER[N_HISTORY]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                             1, 1, 1, 1, 1, 1, 1, 1, 1};
constexpr int HISTORY_WEIGHT_REFUTATION_MASTER[N_HISTORY] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                             0, 0, 0, 0, 0, 0, 0, 0, 0};

constexpr int HISTORY_SCALE_CAPTURE_MAIN_MASTER[N_HISTORY]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 16,
                                                               1, 1, 1, 1, 1, 1, 1, 1, 1};
constexpr int HISTORY_WEIGHT_CAPTURE_MAIN_MASTER[N_HISTORY] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                                               0, 0, 0, 0, 0, 0, 0, 0, 0};


constexpr int HISTORY_SCALE_START[N_HISTORY]  = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                                 1, 1, 1, 1, 1, 1, 1, 1, 1};
constexpr int HISTORY_WEIGHT_START[N_HISTORY] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, -1};

void init_stats(bool onlyD = false);


//---------------
// search + measurement

constexpr int  HISTORY_BUCKETS     = 10000;
constexpr bool USE_DEPTH_WEIGHT    = true;
constexpr bool USE_ONLY_RANK       = true;
constexpr bool USE_FULL_QUIET_SORT = false;
constexpr bool MEASURE_BIAS        = true;

constexpr bool STATS_REFUTATION           = false;
constexpr bool STATS_QUIETS               = true;
constexpr bool STATS_QUIET_EVASION_MAIN   = false;
constexpr bool STATS_QUIET_EVASION_QS     = false;
constexpr bool STATS_CAPTURE_EVASION_MAIN = false;
constexpr bool STATS_CAPTURE_EVASION_QS   = false;
constexpr bool STATS_CAPTURE_MAIN         = false;

static_assert(!(STATS_REFUTATION
                && (STATS_QUIETS || STATS_QUIET_EVASION_MAIN || STATS_QUIET_EVASION_QS
                    || STATS_CAPTURE_EVASION_MAIN || STATS_CAPTURE_EVASION_QS
                    || STATS_CAPTURE_MAIN)));
static_assert(!(STATS_QUIETS
                && (STATS_REFUTATION || STATS_QUIET_EVASION_MAIN || STATS_QUIET_EVASION_QS
                    || STATS_CAPTURE_EVASION_MAIN || STATS_CAPTURE_EVASION_QS
                    || STATS_CAPTURE_MAIN)));
static_assert(!(STATS_QUIET_EVASION_MAIN
                && (STATS_REFUTATION || STATS_QUIETS  //|| STATS_CAPTURE_EVASION_MAIN
                    || STATS_CAPTURE_EVASION_QS || STATS_CAPTURE_MAIN)));
static_assert(!(STATS_QUIET_EVASION_QS
                && (STATS_REFUTATION || STATS_QUIETS  //|| STATS_CAPTURE_EVASION_MAIN
                    || STATS_CAPTURE_EVASION_QS || STATS_CAPTURE_MAIN)));
static_assert(!(STATS_CAPTURE_EVASION_MAIN
                && (STATS_REFUTATION || STATS_QUIETS  //|| STATS_QUIET_EVASION_MAIN
                    || STATS_QUIET_EVASION_QS || STATS_CAPTURE_MAIN)));
static_assert(!(STATS_CAPTURE_EVASION_QS
                && (STATS_REFUTATION || STATS_QUIETS  //|| STATS_QUIET_EVASION_MAIN
                    || STATS_QUIET_EVASION_QS || STATS_CAPTURE_MAIN)));
static_assert(!(STATS_CAPTURE_MAIN
                && (STATS_REFUTATION || STATS_QUIETS || STATS_QUIET_EVASION_MAIN
                    || STATS_QUIET_EVASION_QS || STATS_CAPTURE_EVASION_MAIN
                    || STATS_CAPTURE_EVASION_QS)));

//---------------
// uci stats command

constexpr std::tuple<int, int, const char*> STATS_STEPS[] = {
  /*
  // refutations
  {0, 1, "0"   },
  {1, 1, "1"   },
  */
  // quiets
  //{-2, 1, "-2"   },
  {-1, 1, "-1"   },
  {-1, 2, "-0.5" },
  {-1, 4, "-0.25"},
  {0,  1, "0"    },
  {1,  4, "0.25" },
  {1,  2, "0.5"  },
  {1,  1, "1"    },
 //{2, 1, "2"   },
  /*
  // captures
  {-1, 16, "-1/16"},
  {-1, 32, "-1/32"},
  {-1, 64, "-1/64"},
  {0,  1,  "0"    },
  {1,  64, "1/64" },
  {1,  32, "1/32" },
  {1,  16, "1/16" },
  {2,  16, "1/8"  },
  */
};

constexpr std::tuple<int, const char*> STATS_PARAMS[] = {
  //{HISTORY_CAPTURE, "capture"},
  //{HISTORY_MAIN, "main"},
  //{HISTORY_PAWN, "pawn"},
  //{HISTORY_INCHECK, "incheck"},
  //{HISTORY_CMH0, "cmh0"},
  //{HISTORY_CMH1, "cmh1"},
  //{HISTORY_CMH2, "cmh2"},
  //{HISTORY_CMH3, "cmh3"},
  //{HISTORY_CMH4, "cmh4"},
  //{HISTORY_CMH5, "cmh5"},
  //{HISTORY_GIVES_CHECK, "givesCheck"},
  {HISTORY_ESCAPE,   "escape" },
  {HISTORY_EN_PRISE, "enprise"},
 //{HISTORY_CMH0_POS, "cmh0_pos"},
  //{HISTORY_CMH0_NEG, "cmh0_neg"},
  //{HISTORY_REF_ORDER, "k1k2cm"},
};

inline int getBucket(int V) {
    return std::clamp(int(int64_t(V - Dmin) * HISTORY_BUCKETS / (Dmax - Dmin)), 0,
                      HISTORY_BUCKETS - 1);
}

inline int getWeight(int depth) {
    return USE_DEPTH_WEIGHT
            && !((STATS_QUIET_EVASION_QS && !STATS_QUIET_EVASION_MAIN)
                 || (STATS_CAPTURE_EVASION_QS && !STATS_CAPTURE_EVASION_MAIN))
           ? depth
               * (1
                  + ((STATS_QUIET_EVASION_MAIN && STATS_QUIET_EVASION_QS)
                     || (STATS_CAPTURE_EVASION_MAIN && STATS_CAPTURE_EVASION_QS)))
           : 1;
}

}
#endif
