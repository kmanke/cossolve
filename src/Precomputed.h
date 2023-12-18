/* Here we define precomputed functions and linearizations.
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

#ifndef COSSOLVE_PRECOMPUTED_H
#define COSSOLVE_PRECOMPUTED_H

#include "CossolveTypes.h"

namespace cossolve {

// BCH approximations
// This function computes the Jacobian of the BCH formula to the specified order with respect
// to the specified argument (0 for the first argument, 1 for the second)
template <unsigned order, unsigned argument>
void bchJacobian(Eigen::Ref<const TwistType> x0, Eigen::Ref<const TwistType> y0,
		 SingleMatrixType& J);

// Strain jacobians
void computeStrainJacobian(Eigen::Ref<const TwistType> strain0, Eigen::Ref<const TwistType> freeStrain0,
			   const SingleMatrixType& stiffness, SingleMatrixType& Jf);
void computeStrainInitialTerm(Eigen::Ref<const TwistType> strain0, Eigen::Ref<const TwistType> freeStrain0,
			      const SingleMatrixType& stiffness, Eigen::Ref<TwistType> res);

} // namespace cossolve

#endif // COSSOLVE_PRECOMPUTED_H
