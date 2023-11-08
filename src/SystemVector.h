/* This class represents the main system vector for the solver
 * with its subvectors.
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

#ifndef COSSOLVE_SYSTEM_VECTOR_H
#define COSSOLVE_SYSTEM_VECTOR_H

#include "config.h"
#include "CossolveTypes.h"
#include "ConstraintList.h"
#include "ForceList.h"
#include "SolverParameters.h"

namespace cossolve {

class SystemVector
{
public:
    SystemVector(const SolverParameters& params, const ForceList& appliedForces,
		 const ConstraintList& fixedConstraints);
    ~SystemVector() { }

    Eigen::Ref<const VectorType> getStrain() const;
    Eigen::Ref<const VectorType> getFixedConstraintForces() const;

    // Initializes constants which are used in applied force calculations
    void initAppliedForceConstants(const MatrixType& Einv, const ScalarType fStar0);
    
    // Updates the fixed constraint in the RHS
    void updateFixedConstraints();

    // Updates the applied forces in the RHS
    void updateAppliedForces(const MatrixType& EinvKD, const MatrixType& AfK,
			     const VectorType& fStar, const CoordList& gBody);
    
private:
    enum class SubVector
    {
	lhs,
	strain,
	fixedConstraintForce,
	rhs,
	appliedForce,
	fixedConstraintTwist
    };
    // The full system vector
    // We split into the left hand side (inputs / solver outputs) and rhs (outputs)
    VectorType lhs;
    VectorType rhs;

    // References to some shared structures
    const SolverParameters& params;
    const ForceList& appliedForces;
    const ConstraintList& fixedConstraints;

    // We store some vectors here to avoid recomputing
    VectorType EinvfStar0;

    // Utilities for indexing the subvectors
    template <SubVector sub>
    int firstRow();
    template <SubVector sub>
    int lastRow();
    template <SubVector sub>
    int nRows();
    template <SubVector sub>
    int rowIndex(int index, int stride = 1);
    template <SubVector sub>
    Eigen::Ref<VectorType> subRef(int index, int length));
};

// Definition of templated functions
template <SystemVector::SubVector sub>
int SystemVector::firstRow()
{
    if constexpr (sub == SubVector::lhs ||
		  sub == SubVector::rhs)
    {
	return 0;
    }
    if constexpr (sub == SubVector::strain ||
		  sub == SubVector::appliedForce)
    {
	return 0;
    }
    else if constexpr (sub == SubVector::fixedConstraintForce ||
		       sub == SubVector::fixedConstraintTwist)
    {
	return lastRow<SubVector::strain>() + 1;
    }
    else
    {
	return 0;
    }
}
    
template <SystemVector::SubVector sub>
int SystemVector::lastRow()
{
    return firstRow<sub>() + nRows<sub>() - 1;
}
    
template <SystemVector::SubVector sub>
int SystemVector::nRows()
{
    if constexpr (sub == SubVector::lhs ||
		  sub == SubVector::rhs)
    {
	return nRows<SubVector::strain>() + nRows<SubVector::fixedConstraintForce>();
    }
    if constexpr (sub == SubVector::strain ||
		  sub == SubVector::appliedForce)
    {
	return params.nNodes() * twistLength;
    }
    else if constexpr (sub == SubVector::fixedConstraintForce ||
		       sub == SubVector::fixedConstraintTwist)
    {
	return fixedConstraints.size() * twistLength;
    }
    else
    {
	return 0;
    }
}

template <SystemVector::SubVector sub>
int SystemVector::rowIndex(int index, int stride)
{
    return (index >= 0) ? (firstRow<sub>() + index*stride) : (lastRow<sub>() + 1 + index*stride);
}

template <SystemVector::SubVector sub>
Eigen::Ref<VectorType> SystemVector::subRef(int index, int length)
{
    if constexpr (sub == SubVector::strain ||
		  sub == SubVector::fixedConstraintForce)
    {
	return lhs.block(rowIndex<sub>(index, length), 0, length, 1);
    }
    else if constexpr (sub == SubVector::appliedForce ||
		       sub == SubVector::fixedConstraintTwist)
    {
	return rhs.block(rowIndex<sub>(index, length), 0, length, 1);
    }
    else
    {
	// Unreachable
	return lhs;
    }
}

} // namespace cossolve

#endif // COSSOLVE_SYSTEM_VECTOR_H
