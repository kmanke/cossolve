/* This header contains useful class aliases.
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

#ifndef COSSOLVE_TYPES_H
#define COSSOLVE_TYPES_H

#include <eigen3/Eigen/Eigen>

namespace cossolve {

constexpr int twistLength = 6;
    
using ScalarType = double; // This defines the data type of all numerical quantities used by the solver
using VectorType = Eigen::VectorX<ScalarType>;
using TwistType = Eigen::Vector<ScalarType, twistLength>;
using MatrixType = Eigen::SparseMatrix<ScalarType>;
using SingleMatrixType = Eigen::Matrix<ScalarType, 6, 6>;
using CoordType = Eigen::Matrix<ScalarType, 4, 4>;
using PointType = Eigen::Vector<ScalarType, 3>;

using TripletList = std::vector<Eigen::Triplet<ScalarType>>;
using PairList = std::vector<std::pair<int, ScalarType>>;

using CoordList = std::vector<CoordType>;

} // namespace cossolve

#endif // COSSOLVE_TYPES_H
