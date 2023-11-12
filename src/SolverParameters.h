/* This struct holds all of the globally accessible solver parameters.
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

#ifndef COSSOLVE_SOLVER_PARAMETERS_H
#define COSSOLVE_SOLVER_PARAMETERS_H

#include "config.h"
#include "CossolveTypes.h"

namespace cossolve {

struct SolverParameters
{
public:
    SolverParameters(ScalarType tensileModulus, ScalarType poissonRatio, ScalarType shearModulus,
		     ScalarType momentX, ScalarType momentY, ScalarType momentZ, ScalarType area,
		     ScalarType length, ScalarType linearDensity, int nNodes)
	: tensileModulus(tensileModulus), poissonRatio(poissonRatio),
	  shearModulus(shearModulus == 0 ? tensileModulus / (2*(1 + poissonRatio)) : shearModulus),
	  momentX(momentX), momentY(momentY), momentZ(momentZ), area(area), length(length),
	linearDensity(linearDensity), nNodes(nNodes), nSegments(nNodes - 1),
	ds(length / nSegments)
    {
    }
    // Physical values
    ScalarType tensileModulus;
    ScalarType poissonRatio;
    ScalarType shearModulus;
    ScalarType momentX;
    ScalarType momentY;
    ScalarType momentZ;
    ScalarType area;
    ScalarType length;
    ScalarType linearDensity;

    // Discretization parameters
    int nNodes; // The number of nodes to place on the rod
    int nSegments; // The number of segments along the rod
    int nFixedConstraints; // The number of fixed constraints
    ScalarType ds; // The change in arclength parameter per segment
};
    
} // namespace cossolve

#endif // COSSOLVE_SOLVER_PARAMETERS_H

