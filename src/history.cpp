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

#include "history.h"

namespace Stockfish {

namespace {

constexpr int UPDATE_POWER_SIZE = 30000;  // Have to be >= then the maximum history table divisior D
constexpr int UPDATE_POWER_SCALE                = 1024;
constexpr int UPDATE_POWER_EXPONENT_NOMINATOR   = 255;
constexpr int UPDATE_POWER_EXPONENT_DENOMINATOR = 256;

int UpdatePower[UPDATE_POWER_SIZE + 1];

}

namespace History {

int getUpdate(int entry, int bonus, int D) {

    assert(D <= UPDATE_POWER_SIZE);

    // Make sure that bonus is in range [-D / 2, D / 2]
    int clampedBonus = std::clamp(bonus, -D / 2, D / 2);
    //assert(std::abs(entry) >= 0);
    //assert(std::abs(entry) <= UPDATE_POWER_SIZE);
    //assert(std::abs(clampedBonus) >= 0);
    //assert(std::abs(clampedBonus) <= UPDATE_POWER_SIZE);
    entry += (clampedBonus - entry * std::abs(clampedBonus) / D) * UpdatePower[std::abs(entry)]
         * UPDATE_POWER_EXPONENT_NOMINATOR
         / (UPDATE_POWER_EXPONENT_DENOMINATOR * UpdatePower[std::abs(clampedBonus)]);
    return std::clamp(entry, -D, D);
}

void init() {

    constexpr double A =
      double(UPDATE_POWER_EXPONENT_NOMINATOR) / UPDATE_POWER_EXPONENT_DENOMINATOR;
    constexpr double E = (A - 1) / A;

    UpdatePower[0] = 1;
    for (int i = 1; i <= UPDATE_POWER_SIZE; i++)
        UpdatePower[i] = std::max(int(UPDATE_POWER_SCALE * std::pow(i, E)), 1);
    //for (int i = 0; i <= UPDATE_POWER_SIZE; i++)
    //    std::cerr << i << " " << UpdatePower[i] << std::endl;
    //std::exit(1);
}

}  // namespace History

}  // namespace Stockfish
