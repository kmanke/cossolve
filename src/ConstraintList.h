/* This class represents a list of constraints which is shared
 * between various solver elements.
 *
 * Copyright (C) 2023 Kyle Manke
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef COSSOLVE_CONSTRAINT_LIST_H
#define COSSOLVE_CONSTRAINT_LIST_H

#include "CossolveTypes.h"
#include "config.h"

#include <vector>

namespace cossolve {

struct ConstraintListEntry
{
    int node;
    CoordType g;
    bool active;
};

using ConstraintList = std::vector<ConstraintListEntry>;
    
} // namespace cossolve
    
#endif // COSSOLVE_CONSTRAINT_LIST_H
