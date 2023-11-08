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
#include "SolverParameters.h"

#include <eigen3/Eigen/Eigen>

namespace cossolve {

class SystemVector
{
public:
    SystemVector(const SolverParameters& params, const ConstraintList& fixedConstraints);
    ~SystemVector() { }

    Eigen::Ref<const VectorType> getStrain() const;
    Eigen::Ref<const VectorType> getFixedConstraintForces() const;
    
private:
    // The full system vector
    // We split into the left hand side (inputs / solver outputs) and rhs (outputs)
    VectorType lhs;
    VectorType rhs;

    // References to some shared structures
    const SolverParameters& params;
    const ConstraintList& fixedConstraints;
};

} // namespace cossolve

#endif // COSSOLVE_SYSTEM_VECTOR_H
