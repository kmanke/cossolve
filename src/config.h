/* This header file is used to convert cmake build options to
 * flags and variables that can be used throughout the project.
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

#ifndef COSSOLVE_CONFIG_H
#define COSSOLVE_CONFIG_H

#ifdef USE_MKL
#include <eigen3/Eigen/PardisoSupport>
#else
#include <eigen3/Eigen/UmfPackSupport>
#endif

namespace cossolve {

#ifdef USE_MKL
template <typename T>
using SparseLUSolver = Eigen::PardisoLU<T>;
#else
template <typename T>
using SparseLUSolver = Eigen::UmfPackLU<T>;
#endif

} // namespace cossolve

#endif
