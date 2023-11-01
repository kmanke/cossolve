/* Functions for some common operations on SE(3) / SO(3) and their associated Lie algebras.
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

#ifndef COSSOLVE_SE3_H
#define COSSOLVE_SE3_H

#include <eigen3/Eigen/Eigen>
#include <cmath>

namespace cossolve
{

// Hat operator from R3 to so(3).
// Inputs:
//   x - vector in R3 to transform.
//   xhat - 3x3 matrix where the result will be stored.
void so3hat(Eigen::Ref<const Eigen::Matrix<double, 3, 1>> x, Eigen::Ref<Eigen::Matrix<double, 3, 3>> xhat);

// Overload of so3hat for sparse matrices.
void so3hat(Eigen::Ref<const Eigen::Matrix<double, 3, 1>> x, Eigen::Ref<Eigen::SparseMatrix<double>> xhat);

// Hat operator from R6 to se(3).
// Inputs:
//   x - vector in R6 to transform.
//   xhat - 4x4 matrix where the result will be stored.
void se3hat(Eigen::Ref<const Eigen::Matrix<double, 6, 1>> x, Eigen::Ref<Eigen::Matrix<double, 4, 4>> xhat);

// Unhat operator from so(3) to R3.
// Inputs:
//   xhat - matrix in so(3) to transform.
//   x - vector in R3 where the result will be stored.
void so3unhat(Eigen::Ref<const Eigen::Matrix<double, 3, 3>> xhat, Eigen::Ref<Eigen::Matrix<double, 3, 1>> x);

// Unhat operator from se(3) to R6.
// Inputs:
//   xhat - matrix in se(3) to transform.
//   x - vector in R6 where the result will be stored.
void se3unhat(Eigen::Ref<const Eigen::Matrix<double, 4, 4>> xhat, Eigen::Ref<Eigen::Matrix<double, 6, 1>> x);

// Adjoint operator for SE(3).
// Inputs:
//   g - matrix in SE(3) to take the adjoint.
//   Adg - resulting 6x6 Adjoint matrix.
void Adjoint(Eigen::Ref<const Eigen::Matrix<double, 4, 4>> g, Eigen::Ref<Eigen::Matrix<double, 6, 6>> Adg);

// adjoint operator for R6.
// Inputs:
//   x - twist in R6.
//   adx - resulting 6x6 adjoint matrix.
void adjoint(Eigen::Ref<const Eigen::Matrix<double, 6, 1>> x, Eigen::Ref<Eigen::Matrix<double, 6, 6>> adx);

// Overload of adjoint for sparse matrices.
void adjoint(Eigen::Ref<const Eigen::Matrix<double, 6, 1>> x, Eigen::Ref<Eigen::SparseMatrix<double>> adx);

} // namespace cossolve

#endif // COSSOLVE_SE3_H
