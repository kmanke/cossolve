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

class SolverParameters
{
public:
    SolverParameters(ScalarType tensileModulus, ScalarType poissonRatio, ScalarType shearModulus,
		     ScalarType momentX, ScalarType momentY, ScalarType momentZ, ScalarType area,
		     ScalarType length, ScalarType linearDensity, int nNodes)
	: tensileModulusVal(tensileModulus), poissonRatioVal(poissonRatio),
	  shearModulusVal(shearModulus == 0 ? tensileModulus / (2*(1 + poissonRatio)) : shearModulus),
	  momentXVal(momentX), momentYVal(momentY), momentZVal(momentZ), areaVal(area), lengthVal(length),
	linearDensityVal(linearDensity), nNodesVal(nNodes), nSegmentsVal(nNodesVal - 1),
	dsVal(lengthVal / nSegmentsVal)
    {
    }
    
    // Public const references to values
    ScalarType tensileModulus() const { return tensileModulusVal; }
    ScalarType poissonRatio() const { return poissonRatioVal; }
    ScalarType shearModulus() const { return shearModulusVal; }
    ScalarType momentX() const { return momentXVal; }
    ScalarType momentY() const { return momentYVal; }
    ScalarType momentZ() const { return momentZVal; }
    ScalarType area() const { return areaVal; }
    ScalarType length() const { return lengthVal; }
    ScalarType linearDensity() const { return linearDensityVal; }

    int nNodes() const { return nNodesVal; }
    int nSegments() const { return nSegmentsVal; }
    ScalarType ds() const { return dsVal; }
    
private:
    // Physical values
    ScalarType tensileModulusVal;
    ScalarType poissonRatioVal;
    ScalarType shearModulusVal;
    ScalarType momentXVal;
    ScalarType momentYVal;
    ScalarType momentZVal;
    ScalarType areaVal;
    ScalarType lengthVal;
    ScalarType linearDensityVal;

    // Discretization parameters
    int nNodesVal; // The number of nodes to place on the rod
    int nSegmentsVal; // The number of segments along the rod
    int nFixedConstraintsVal; // The number of fixed constraints
    ScalarType dsVal; // The change in arclength parameter per segment
};
    
} // namespace cossolve

#endif // COSSOLVE_SOLVER_PARAMETERS_H

