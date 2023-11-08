/* Implementation of SystemVector
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

#include "SystemVector.h"
#include "LieGroupOps.h"

#include <eigen3/unsupported/Eigen/MatrixFunctions>

namespace cossolve {

SystemVector::SystemVector(const SolverParameters& params, const ForceList& appliedForces,
			   const ConstraintList& fixedConstraints)
    : params(params), appliedForces(appliedForces), fixedConstraints(fixedConstraints)
{
    // Allocate our vectors
    lhs = VectorType(nRows<SubVector::lhs>());
    rhs = VectorType(nRows<SubVector::rhs>());

    return;
}

void SystemVector::initAppliedForceConstants(const MatrixType& Einv, const ScalarType fStar0)
{
    EinvfStar0 = VectorType::Ones(nRows<SubVector::appliedForce>());
    EinvfStar0 *= (2.0 / params.ds());
    EinvfStar0 = Einv * EinvfStar0;

    return;
}

void SystemVector::updateAppliedForces(const MatrixType& EinvKD, const MatrixType& AfK,
				       const VectorType& fStar, const CoordList& gBody)
{
    Eigen::Ref<VectorType> forces = subRef<SubVector::appliedForce>(0, nRows<SubVector::appliedForce>());
    forces = -AfK*fStar + EinvfStar0;

    // Iterate through the force list and add all of the forces
    CoordType g;
    SingleMatrixType Adg;
    
    // Apply all external forces
    for (auto it = appliedForces.cbegin(); it != appliedForces.cend(); ++it)
    {	
	if (it->bodyFrame) // Force can be applied directly with no transformation
	{
	    subRef<SubVector::appliedForce>(it->node, twistLength) += it->force;
	}
	else // Need to rotate from the spatial frame into the body frame
	{
	    // Construct a coordinate transformation which only rotates
	    g.setIdentity();
	    g.block<3, 3>(0, 0) = gBody[it->node].block<3, 3>(0, 0);
	    Adg.setZero();
	    Adjoint(g, Adg);
	    subRef<SubVector::appliedForce>(it->node, twistLength) +=
		Adg.transpose() * it->force;
	}
    }
    return; 

}

void SystemVector::updateFixedConstraints()
{
    Eigen::Ref<VectorType> constraintTwists = subRef<SubVector::fixedConstraintTwist>(0, twistLength);
    constraintTwists.setZero();

    // Iterate through the constraints, only filling in the active ones
    int index = 0;
    TwistType twist;
    CoordType twistHat;
    for (auto it = fixedConstraints.cbegin(); it != fixedConstraints.cend(); ++it)
    {
	if (it->active)
	{
	    twist.setZero();
	    twistHat = it->g.log();
	    se3unhat(twistHat, twist);
	    subRef<SubVector::fixedConstraintTwist>(index, twistLength);
	}
	index++;
    }
    return;
}

} // namespace cossolve
