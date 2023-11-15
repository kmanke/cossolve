/* Here we define functions to calculate the Jacobians for estimating
 * twist as a function of phi.
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

#ifndef COSSOLVE_TWIST_JACOBIAN_H
#define COSSOLVE_TWIST_JACOBIAN_H

#include "CossolveTypes.h"

namespace cossolve {

// Twist jacobians
void computePhiJacobian(Eigen::Ref<const TwistType> ksi0, Eigen::Ref<const TwistType> phi0,
			SingleMatrixType& Jphi);
void computeTwistJacobian(Eigen::Ref<const TwistType> ksi0, Eigen::Ref<const TwistType> phi0,
			  SingleMatrixType& Jksi);
void computeInitialTerm(Eigen::Ref<const TwistType> ksi0, Eigen::Ref<const TwistType> phi0, Eigen::Ref<TwistType> res);

// Strain jacobians
void computeStrainJacobian(Eigen::Ref<const TwistType> strain0, Eigen::Ref<const TwistType> freeStrain0,
			   const SingleMatrixType& stiffness, SingleMatrixType& Jf);
void computeStrainInitialTerm(Eigen::Ref<const TwistType> strain0, Eigen::Ref<const TwistType> freeStrain0,
			      const SingleMatrixType& stiffness, Eigen::Ref<TwistType> res);

} // namespace cossolve

#endif // COSSOLVE_TWIST_JACOBIAN_H
