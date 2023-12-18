/* Implementation of TwistJacobian functions.
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

#include "Precomputed.h"

namespace cossolve {

template <>
void bchJacobian<3, 0>(Eigen::Ref<const TwistType> _x0, Eigen::Ref<const TwistType> _y0,
		       SingleMatrixType& J)
{
    ScalarType x0 = _x0(0), x1 = _x0(1), x2 = _x0(2), x3 = _x0(3), x4 = _x0(4), x5 = _x0(5);
    ScalarType y0 = _y0(0), y1 = _y0(1), y2 = _y0(2), y3 = _y0(3), y4 = _y0(4), y5 = _y0(5);

    J(0, 0) = x4*y4 + x5*y5 - y4*y4 - y5*y5 + 12;
    J(0, 1) = x3*y4 - 2*x4*y3 + y3*y4 + 6*y5;
    J(0, 2) = x3*y5 - 2*x5*y3 + y3*y5 - 6*y4;
    J(0, 3) = x1*y4 + x2*y5 + x4*y1 + x5*y2 - 2*y1*y4 - 2*y2*y5;
    J(0, 4) = x0*y4 - 2*x1*y3 + x3*y1 - 2*x4*y0 + y0*y4 + y1*y3 + 6*y2;
    J(0, 5) = x0*y5 - 2*x2*y3 + x3*y2 - 2*x5*y0 + y0*y5 - 6*y1 + y2*y3;
    
    J(1, 0) = -2*x3*y4 + x4*y3 + y3*y4 - 6*y5;
    J(1, 1) = x3*y3 + x5*y5 - y3*y3 - y5*y5 + 12;
    J(1, 2) = x4*y5 - 2*x5*y4 + 6*y3 + y4*y5;
    J(1, 3) = -2*x0*y4 + x1*y3 - 2*x3*y1 + x4*y0 + y0*y4 + y1*y3 - 6*y2;
    J(1, 4) = x0*y3 + x2*y5 + x3*y0 + x5*y2 - 2*y0*y3 - 2*y2*y5;
    J(1, 5) = x1*y5 - 2*x2*y4 + x4*y2 - 2*x5*y1 + 6*y0 + y1*y5 + y2*y4;
    
    J(2, 0) = -2*x3*y5 + x5*y3 + y3*y5 + 6*y4;
    J(2, 1) = -2*x4*y5 + x5*y4 - 6*y3 + y4*y5;
    J(2, 2) = x3*y3 + x4*y4 - y3*y3 - y4*y4 + 12;
    J(2, 3) = -2*x0*y5 + x2*y3 - 2*x3*y2 + x5*y0 + y0*y5 + 6*y1 + y2*y3;
    J(2, 4) = -2*x1*y5 + x2*y4 - 2*x4*y2 + x5*y1 - 6*y0 + y1*y5 + y2*y4;
    J(2, 5) = x0*y3 + x1*y4 + x3*y0 + x4*y1 - 2*y0*y3 - 2*y1*y4;
    
    J(3, 0) = 0;
    J(3, 1) = 0;
    J(3, 2) = 0;
    J(3, 3) = x4*y4 + x5*y5 - y4*y4 - y5*y5 + 12;
    J(3, 4) = x3*y4 - 2*x4*y3 + y3*y4 + 6*y5;
    J(3, 5) = x3*y5 - 2*x5*y3 + y3*y5 - 6*y4;
    
    J(4, 0) = 0;
    J(4, 1) = 0;
    J(4, 2) = 0;
    J(4, 3) = -2*x3*y4 + x4*y3 + y3*y4 - 6*y5;
    J(4, 4) = x3*y3 + x5*y5 - y3*y3 - y5*y5 + 12;
    J(4, 5) = x4*y5 - 2*x5*y4 + 6*y3 + y4*y5;

    J(5, 0) = 0;
    J(5, 1) = 0;
    J(5, 2) = 0;
    J(5, 3) = -2*x3*y5 + x5*y3 + y3*y5 + 6*y4;
    J(5, 4) = -2*x4*y5 + x5*y4 - 6*y3 + y4*y5;
    J(5, 5) = x3*y3 + x4*y4 - y3*y3 - y4*y4 + 12;

    J /= 12;
    return;
}

template <>
void bchJacobian<3, 1>(Eigen::Ref<const TwistType> _x0, Eigen::Ref<const TwistType> _y0,
		       SingleMatrixType& J)
{
    ScalarType x0 = _x0(0), x1 = _x0(1), x2 = _x0(2), x3 = _x0(3), x4 = _x0(4), x5 = _x0(5);
    ScalarType y0 = _y0(0), y1 = _y0(1), y2 = _y0(2), y3 = _y0(3), y4 = _y0(4), y5 = _y0(5);

    J(0, 0) = -x4*x4 + x4*y4 - x5*x5 + x5*y5 + 12;
    J(0, 1) = x3*x4 - 2*x3*y4 + x4*y3 - 6*x5;
    J(0, 2) = x3*x5 - 2*x3*y5 + 6*x4 + x5*y3;
    J(0, 3) = -2*x1*x4 + x1*y4 - 2*x2*x5 + x2*y5 + x4*y1 + x5*y2;
    J(0, 4) = x0*x4 - 2*x0*y4 + x1*x3 + x1*y3 - 6*x2 - 2*x3*y1 + x4*y0;
    J(0, 5) = x0*x5 - 2*x0*y5 + 6*x1 + x2*x3 + x2*y3 - 2*x3*y2 + x5*y0;
    
    J(1, 0) = x3*x4 + x3*y4 - 2*x4*y3 + 6*x5;
    J(1, 1) = -x3*x3 + x3*y3 - x5*x5 + x5*y5 + 12;
    J(1, 2) = -6*x3 + x4*x5 - 2*x4*y5 + x5*y4;
    J(1, 3) = x0*x4 + x0*y4 + x1*x3 - 2*x1*y3 + 6*x2 + x3*y1 - 2*x4*y0;
    J(1, 4) = -2*x0*x3 + x0*y3 - 2*x2*x5 + x2*y5 + x3*y0 + x5*y2;
    J(1, 5) = -6*x0 + x1*x5 - 2*x1*y5 + x2*x4 + x2*y4 - 2*x4*y2 + x5*y1;

    J(2, 0) = x3*x5 + x3*y5 - 6*x4 - 2*x5*y3;
    J(2, 1) = 6*x3 + x4*x5 + x4*y5 - 2*x5*y4;
    J(2, 2) = -x3*x3 + x3*y3 - x4*x4 + x4*y4 + 12;
    J(2, 3) = x0*x5 + x0*y5 - 6*x1 + x2*x3 - 2*x2*y3 + x3*y2 - 2*x5*y0;
    J(2, 4) = 6*x0 + x1*x5 + x1*y5 + x2*x4 - 2*x2*y4 + x4*y2 - 2*x5*y1;
    J(2, 5) = -2*x0*x3 + x0*y3 - 2*x1*x4 + x1*y4 + x3*y0 + x4*y1;

    J(3, 0) = 0;
    J(3, 1) = 0;
    J(3, 2) = 0;
    J(3, 3) = -x4*x4 + x4*y4 - x5*x5 + x5*y5 + 12;
    J(3, 4) = x3*x4 - 2*x3*y4 + x4*y3 - 6*x5;
    J(3, 5) = x3*x5 - 2*x3*y5 + 6*x4 + x5*y3;

    J(4, 0) = 0;
    J(4, 1) = 0;
    J(4, 2) = 0;
    J(4, 3) = x3*x4 + x3*y4 - 2*x4*y3 + 6*x5;
    J(4, 4) = -x3*x3 + x3*y3 - x5*x5 + x5*y5 + 12;
    J(4, 5) = -6*x3 + x4*x5 - 2*x4*y5 + x5*y4;

    J(5, 0) = 0;
    J(5, 1) = 0;
    J(5, 2) = 0;
    J(5, 3) = x3*x5 + x3*y5 - 6*x4 - 2*x5*y3;
    J(5, 4) = 6*x3 + x4*x5 + x4*y5 - 2*x5*y4;
    J(5, 5) = -x3*x3 + x3*y3 - x4*x4 + x4*y4 + 12;

    J /= 12;
    return;
}

template <>
void bchJacobian<2, 0>(Eigen::Ref<const TwistType> _x0, Eigen::Ref<const TwistType> _y0,
		       SingleMatrixType& J)
{
    ScalarType y0 = _y0(0), y1 = _y0(1), y2 = _y0(2), y3 = _y0(3), y4 = _y0(4), y5 = _y0(5);

    J(0, 0) = 2;
    J(0, 1) = y5;
    J(0, 2) = -y4;
    J(0, 3) = 0;
    J(0, 4) = y2;
    J(0, 5) = -y1;
    
    J(1, 0) = -y5;
    J(1, 1) = 2;
    J(1, 2) = y3;
    J(1, 3) = -y2;
    J(1, 4) = 0;
    J(1, 5) = y0;
    
    J(2, 0) = y4;
    J(2, 1) = -y3;
    J(2, 2) = 2;
    J(2, 3) = y1;
    J(2, 4) = -y0;
    J(2, 5) = 0;

    J(3, 0) = 0;
    J(3, 1) = 0;
    J(3, 2) = 0;
    J(3, 3) = 2;
    J(3, 4) = y5;
    J(3, 5) = -y4;
    
    J(4, 0) = 0;
    J(4, 1) = 0;
    J(4, 2) = 0;
    J(4, 3) = -y5;
    J(4, 4) = 2;
    J(4, 5) = y3;
    
    J(5, 0) = 0;
    J(5, 1) = 0;
    J(5, 2) = 0;
    J(5, 3) = y4;
    J(5, 4) = -y3;
    J(5, 5) = 2;

    J /= 2;
    return;
}

template <>
void bchJacobian<2, 1>(Eigen::Ref<const TwistType> _x0, Eigen::Ref<const TwistType> _y0,
		       SingleMatrixType& J)
{
    ScalarType x0 = _x0(0), x1 = _x0(1), x2 = _x0(2), x3 = _x0(3), x4 = _x0(4), x5 = _x0(5);

    J(0, 0) = 2;
    J(0, 1) = -x5;
    J(0, 2) = x4;
    J(0, 3) = 0;
    J(0, 4) = -x2;
    J(0, 5) = x1;

    J(1, 0) = x5;
    J(1, 1) = 2;
    J(1, 2) = -x3;
    J(1, 3) = x2;
    J(1, 4) = 0;
    J(1, 5) = -x0;
    
    J(2, 0) = -x4;
    J(2, 1) = x3;
    J(2, 2) = 2;
    J(2, 3) = -x1;
    J(2, 4) = x0;
    J(2, 5) = 0;
    
    J(3, 0) = 0;
    J(3, 1) = 0;
    J(3, 2) = 0;
    J(3, 3) = 2;
    J(3, 4) = -x5;
    J(3, 5) = x4;
    
    J(4, 0) = 0;
    J(4, 1) = 0;
    J(4, 2) = 0;
    J(4, 3) = x5;
    J(4, 4) = 2;
    J(4, 5) = -x3;
    
    J(5, 0) = 0;
    J(5, 1) = 0;
    J(5, 2) = 0;
    J(5, 3) = -x4;
    J(5, 4) = x3;
    J(5, 5) = 2;

    J /= 2;
    return;
}

template <>
void bchJacobian<1, 0>(Eigen::Ref<const TwistType> _x0, Eigen::Ref<const TwistType> _y0,
		       SingleMatrixType& J)
{
    J.setIdentity();
    return;
};

template <>
void bchJacobian<1, 1>(Eigen::Ref<const TwistType> _x0, Eigen::Ref<const TwistType> _y0,
		       SingleMatrixType& J)
{
    J.setIdentity();
    return;
};

void computeStrainJacobian(Eigen::Ref<const TwistType> strain0, Eigen::Ref<const TwistType> freeStrain0,
			   const SingleMatrixType& stiffness, SingleMatrixType& Jf)
{
    ScalarType f0 = strain0(0), f1 = strain0(1), f2 = strain0(2), f3 = strain0(3), f4 = strain0(4), f5 = strain0(5);
    ScalarType fs0 = freeStrain0(0), fs1 = freeStrain0(1), fs2 = freeStrain0(2), fs3 = freeStrain0(3),
	fs4 = freeStrain0(4), fs5 = freeStrain0(5);
    ScalarType k0 = stiffness(0, 0), k1 = stiffness(1, 1), k2 = stiffness(2, 2), k3 = stiffness(3, 3),
	k4 = stiffness(4, 4), k5 = stiffness(5, 5);
    // dg/df_0
    Jf(0, 0) = 0;
    Jf(0, 1) = f5*k1;
    Jf(0, 2) = -f4*k2;
    Jf(0, 3) = 0;
    Jf(0, 4) = k2*(-f2 + fs2);
    Jf(0, 5) = k1*(f1 - fs1);
    // dg/df_1
    Jf(1, 0) = -f5*k0;
    Jf(1, 1) = 0;
    Jf(1, 2) = f3*k2;
    Jf(1, 3) = k2*(f2 - fs2);
    Jf(1, 4) = 0;
    Jf(1, 5) = k0*(-f0 + fs0);
    // dg/df_2
    Jf(2, 0) = f4*k0;
    Jf(2, 1) = -f3*k1;
    Jf(2, 2) = 0;
    Jf(2, 3) = k1*(-f1 + fs1);
    Jf(2, 4) = k0*(f0 - fs0);
    Jf(2, 5) = 0;
    // dg/df_3
    Jf(3, 0) = 0;
    Jf(3, 1) = f2*k1 - k2*(f2 - fs2);
    Jf(3, 2) = -f1*k2 + k1*(f1 - fs1);
    Jf(3, 3) = 0;
    Jf(3, 4) = f5*k4 - k5*(f5 - fs5);
    Jf(3, 5) = -f4*k5 + k4*(f4 - fs4);
    // dg/df_4
    Jf(4, 0) = -f2*k0 + k2*(f2 - fs2);
    Jf(4, 1) = 0;
    Jf(4, 2) = f0*k2 - k0*(f0 - fs0);
    Jf(4, 3) = -f5*k3 + k5*(f5 - fs5);
    Jf(4, 4) = 0;
    Jf(4, 5) = f3*k5 - k3*(f3 - fs3);
    // dg/df_5
    Jf(5, 0) = f1*k0 - k1*(f1 - fs1);
    Jf(5, 1) = -f0*k1 + k0*(f0 - fs0);
    Jf(5, 2) = 0;
    Jf(5, 3) = f4*k3 - k4*(f4 - fs4);
    Jf(5, 4) = -f3*k4 + k3*(f3 - fs3);
    Jf(5, 5) = 0;
    
    return;
}

void computeStrainInitialTerm(Eigen::Ref<const TwistType> strain0, Eigen::Ref<const TwistType> freeStrain0,
			      const SingleMatrixType& stiffness, Eigen::Ref<TwistType> res)
{
    ScalarType f0 = strain0(0), f1 = strain0(1), f2 = strain0(2), f3 = strain0(3), f4 = strain0(4), f5 = strain0(5);
    ScalarType fs0 = freeStrain0(0), fs1 = freeStrain0(1), fs2 = freeStrain0(2), fs3 = freeStrain0(3),
	fs4 = freeStrain0(4), fs5 = freeStrain0(5);
    ScalarType k0 = stiffness(0, 0), k1 = stiffness(1, 1), k2 = stiffness(2, 2), k3 = stiffness(3, 3),
	k4 = stiffness(4, 4), k5 = stiffness(5, 5);
    
    res(0) = -f4*k2*(f2 - fs2) + f5*k1*(f1 - fs1);
    res(1) = f3*k2*(f2 - fs2) - f5*k0*(f0 - fs0);
    res(2) = -f3*k1*(f1 - fs1) + f4*k0*(f0 - fs0);
    res(3) = -f1*k2*(f2 - fs2) + f2*k1*(f1 - fs1) - f4*k5*(f5 - fs5) + f5*k4*(f4 - fs4);
    res(4) = f0*k2*(f2 - fs2) - f2*k0*(f0 - fs0) + f3*k5*(f5 - fs5) - f5*k3*(f3 - fs3);
    res(5) = -f0*k1*(f1 - fs1) + f1*k0*(f0 - fs0) - f3*k4*(f4 - fs4) + f4*k3*(f3 - fs3);

    return;
}

} // namespace cossolve
