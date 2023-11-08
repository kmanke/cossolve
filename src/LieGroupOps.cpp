/* Implementation of SE(3) / SO(3) functions from se3.h
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

#include "LieGroupOps.h"

#include <iostream>

namespace cossolve
{

void so3hat(Eigen::Ref<const Eigen::Matrix<double, 3, 1>> x, Eigen::Ref<Eigen::Matrix<double, 3, 3>> xhat)
{
    xhat(1, 0) = x(2);
    xhat(0, 1) = -x(2);
    xhat(2, 0) = -x(1);
    xhat(0, 2) = x(1);
    xhat(1, 2) = -x(0);
    xhat(2, 1) = x(0);

    return;
}

void se3hat(Eigen::Ref<const Eigen::Matrix<double, 6, 1>> x, Eigen::Ref<Eigen::Matrix<double, 4, 4>> xhat)
{
    so3hat(x.block<3, 1>(3, 0), xhat.block<3, 3>(0, 0));
    xhat.block<3, 1>(0, 3) = x.block<3, 1>(0, 0);
    
    return;
}

void so3unhat(Eigen::Ref<const Eigen::Matrix<double, 3, 3>> xhat, Eigen::Ref<Eigen::Matrix<double, 3, 1>> x)
{
    x(0) = xhat(2, 1);
    x(1) = xhat(0, 2);
    x(2) = xhat(1, 0);
    
    return;
}

void se3unhat(Eigen::Ref<const Eigen::Matrix<double, 4, 4>> xhat, Eigen::Ref<Eigen::Matrix<double, 6, 1>> x)
{
    so3unhat(xhat.block<3, 3>(0, 0), x.block<3, 1>(3, 0));
    x.block<3, 1>(0, 0) = xhat.block<3, 1>(0, 3);	

    return;
}

void Adjoint(Eigen::Ref<const Eigen::Matrix<double, 4, 4>> g, Eigen::Ref<Eigen::Matrix<double, 6, 6>> Adg)
{
    Eigen::Matrix<double, 3, 3> phat = Eigen::Matrix<double, 3, 3>::Zero();
    so3hat(g.block<3, 1>(0, 3), phat);
    
    Adg.block<3, 3>(0, 0) = g.block<3, 3>(0, 0);
    Adg.block<3, 3>(3, 3) = Adg.block<3, 3>(0, 0);
    Adg.block<3, 3>(0, 3) = phat * g.block<3, 3>(0, 0);

    return;
}

void adjoint(Eigen::Ref<const Eigen::Matrix<double, 6, 1>> x, Eigen::Ref<Eigen::Matrix<double, 6, 6>> adx)
{
    so3hat(x.block<3, 1>(3, 0), adx.block<3, 3>(0, 0));
    so3hat(x.block<3, 1>(0, 0), adx.block<3, 3>(0, 3));
    adx.block<3, 3>(3, 3) = adx.block<3, 3>(0, 0);

    return;
}

} // namespace cossolve
