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

#include "TwistJacobian.h"

namespace cossolve {

void computePhiJacobian(Eigen::Ref<const TwistType> ksi0, Eigen::Ref<const TwistType> phi0,
			SingleMatrixType& Jphi)
{
    ScalarType k0 = ksi0(0), k1 = ksi0(1), k2 = ksi0(2), k3 = ksi0(3), k4 = ksi0(4), k5 = ksi0(5);
    ScalarType p0 = phi0(0), p1 = phi0(1), p2 = phi0(2), p3 = phi0(3), p4 = phi0(4), p5 = phi0(5);
    // df/dp_0
    Jphi(0, 0) = -k4*k4 + k4*p4 - k5*k5 + k5*p5 + 12;
    Jphi(0, 1) = k3*k4 - 2*k3*p4 + k4*p3 - 6*k5;
    Jphi(0, 2) = k3*k5 - 2*k3*p5 + 6*k4 + k5*p3;
    Jphi(0, 3) = -2*k1*k4 + k1*p4 - 2*k2*k5 + k2*p5 + k4*p1 + k5*p2;
    Jphi(0, 4) = k0*k4 - 2*k0*p4 + k1*k3 + k1*p3 - 6*k2 - 2*k3*p1 + k4*p0;
    Jphi(0, 5) = k0*k5 - 2*k0*p5 + 6*k1 + k2*k3 + k2*p3 - 2*k3*p2 + k5*p0;
    // df/dp_1
    Jphi(1, 0) = 3*k4 + k3*p4 - 2*k4*p3 + 6*k5;
    Jphi(1, 1) = -k3*k3 + k3*p3 - k5*k5 + k5*p5 + 12;
    Jphi(1, 2) = -6*k3 + k4*k5 - 2*k4*p5 + k5*p4;
    Jphi(1, 3) = k0*k4 + k0*p4 + k1*k3 - 2*k1*p3 + 6*k2 + k3*p1 - 2*k4*p0;
    Jphi(1, 4) = -2*k0*k3 + k0*p3 - 2*k2*k5 + k2*p5 + k3*p0 + k5*p2;
    Jphi(1, 5) = -6*k0 + k1*k5 - 2*k1*p5 + k2*k4 + k2*p4 - 2*k4*p2 + k5*p1;
    // df/dp_2
    Jphi(2, 0) = k3*k5 + k3*p5 - 6*k4 - 2*k5*p3;
    Jphi(2, 1) = 6*k3 + k4*k5 + k4*p5 - 2*k5*p4;
    Jphi(2, 2) = -k3*k3 + k3*p3 - k4*k4 + k4*p4 + 12;
    Jphi(2, 3) = k0*k5 + k0*p5 - 6*k1 + k2*k3 - 2*k2*p3 + k3*p2 - 2*k5*p0;
    Jphi(2, 4) = 6*k0 + k1*k5 + k1*p5 + k2*k4 - 2*k2*p4 + k4*p2 - 2*k5*p1;
    Jphi(2, 5) = -2*k0*k3 + k0*p3 - 2*k1*k4 + k1*p4 + k3*p0 + k4*p1;
    // df/dp_3
    Jphi(3, 0) = 0;
    Jphi(3, 1) = 0;
    Jphi(3, 2) = 0;
    Jphi(3, 3) = -k4*k4 + k4*p4 - k5*k5 + k5*p5 + 12;
    Jphi(3, 4) = k3*k4 - 2*k3*p4 + k4*p3 - 6*k5;
    Jphi(3, 5) = k3*k5 - 2*k3*p5 + 6*k4 + k5*p3;
    // df/dp_4
    Jphi(4, 0) = 0;
    Jphi(4, 1) = 0;
    Jphi(4, 2) = 0;
    Jphi(4, 3) = k3*k4 + k3*p4 - 2*k4*p3 + 6*k5;
    Jphi(4, 4) = -k3*k3 + k3*p3 - k5*k5 + k5*p5 + 12;
    Jphi(4, 5) = -6*k3 + k4*k5 - 2*k4*p5 + k5*p4;
    // df/dp_5
    Jphi(5, 0) = 0;
    Jphi(5, 1) = 0;
    Jphi(5, 2) = 0;
    Jphi(5, 3) = k3*k5 + k3*p5 - 6*k4 - 2*k5*p3;
    Jphi(5, 4) = 6*k3 + k4*k5 + k4*p5 - 2*k5*p4;
    Jphi(5, 5) = -k3*k3 + k3*p3 - k4*k4 + k4*p4 + 12;

    Jphi /= 12;
    return;
}

void computeTwistJacobian(Eigen::Ref<const TwistType> ksi0, Eigen::Ref<const TwistType> phi0,
			  SingleMatrixType& Jksi)
{
    ScalarType k0 = ksi0(0), k1 = ksi0(1), k2 = ksi0(2), k3 = ksi0(3), k4 = ksi0(4), k5 = ksi0(5);
    ScalarType p0 = phi0(0), p1 = phi0(1), p2 = phi0(2), p3 = phi0(3), p4 = phi0(4), p5 = phi0(5);
    // df/dk_0
    Jksi(0, 0) = -p4*(-k4 + p4) + p5*(k5 - p5);
    Jksi(0, 1) = k3*p4 + p3*(-2*k4 + p4) + 6*p5;
    Jksi(0, 2) = k3*p5 + p3*(-2*k5 + p5) - 6*p4;
    Jksi(0, 3) = -p1*(-k4 + p4) + p2*(k5 - p5) - p4*(-k1 + p1) + p5*(k2 - p2);
    Jksi(0, 4) = k0*p4 + k3*p1 + p0*(-2*k4 + p4) + 6*p2 + p3*(-2*k1 + p1);
    Jksi(0, 5) = k0*p5 + k3*p2 + p0*(-2*k5 + p5) - 6*p1 + p3*(-2*k2 + p2);
    // df/dk_1
    Jksi(1, 0) = k4*p3 + p4*(-2*k3 + p3) - 6*p5;
    Jksi(1, 1) = p3*(k3 - p3) - p5*(-k5 + p5);
    Jksi(1, 2) = k4*p5 + 6*p3 + p4*(-2*k5 + p5);
    Jksi(1, 3) = k1*p3 + k4*p0 + p1*(-2*k3 + p3) - 6*p2 + p4*(-2*k0 + p0);
    Jksi(1, 4) = p0*(k3 - p3) - p2*(-k5 + p5) + p3*(k0 - p0) - p5*(-k2 + p2);
    Jksi(1, 5) = k1*p5 + k4*p2 + 6*p0 + p1*(-2*k5 + p5) + p4*(-2*k2 + p2);
    // df/dk_2
    Jksi(2, 0) = k5*p3 + 6*p4 - p5*(2*k3 - p3);
    Jksi(2, 1) = k5*p4 - 6*p3 - p5*(2*k4 - p4);
    Jksi(2, 2) = p3*(k3 - p3) + p4*(k4 - p4);
    Jksi(2, 3) = k2*p3 + k5*p0 + 6*p1 - p2*(2*k3 - p3) - p5*(2*k0 - p0);
    Jksi(2, 4) = k2*p4 + k5*p1 - 6*p0 - p2*(2*k4 - p4) - p5*(2*k1 - p1);
    Jksi(2, 5) = p0*(k3 - p3) + p1*(k4 - p4) + p3*(k0 - p0) + p4*(k1 - p1);
    // df/dk_3
    Jksi(3, 0) = 0;
    Jksi(3, 1) = 0;
    Jksi(3, 2) = 0;
    Jksi(3, 3) = p4*(k4 - p4) + p5*(k5 - p5);
    Jksi(3, 4) = k3*p4 - p3*(2*k4 - p4) + 6*p5;
    Jksi(3, 5) = k3*p5 - p3*(2*k5 - p5) - 6*p4;
    // df/dk_4
    Jksi(4, 0) = 0;
    Jksi(4, 1) = 0;
    Jksi(4, 2) = 0;
    Jksi(4, 3) = k4*p3 + p4*(-2*k3 + p3) - 6*p5;
    Jksi(4, 4) = p3*(k3 - p3) - p5*(-k5 + p5);
    Jksi(4, 5) = k4*p5 + 6*p3 + p4*(-2*k5 + p5);
    // df/dk_5
    Jksi(5, 0) = 0;
    Jksi(5, 1) = 0;
    Jksi(5, 2) = 0;
    Jksi(5, 3) = k5*p3 + 6*p4 + p5*(-2*k3 + p3);
    Jksi(5, 4) = k5*p4 - 6*p3 + p5*(-2*k4 + p4);
    Jksi(5, 5) = -p3*(-k3 + p3) + p4*(k4 - p4);

    Jksi /= 12;
    return;
}

void computeInitialTerm(Eigen::Ref<const TwistType> ksi0, Eigen::Ref<const TwistType> phi0, Eigen::Ref<TwistType> res)
{
    // Precomputed using numeric solver
    ScalarType k0 = ksi0(0), k1 = ksi0(1), k2 = ksi0(2), k3 = ksi0(3), k4 = ksi0(4), k5 = ksi0(5);
    ScalarType p0 = phi0(0), p1 = phi0(1), p2 = phi0(2), p3 = phi0(3), p4 = phi0(4), p5 = phi0(5);

    res(0) = p0*(-k4*k4 + k4*p4 - k5*k5 + k5*p5 + 12) - p1*(-k3*k4 + k3*p4 + 6*k5) + p2*(k3*k5 - k3*p5 + 6*k4) +
	p3*(-2*k1*k4 + k1*p4 - 2*k2*k5 + k2*p5 + k4*p1 + k5*p2) - p4*(-k0*k4 + k0*p4 - k1*k3 + 6*k2 + k3*p1) +
	p5*(k0*k5 - k0*p5 + 6*k1 + k2*k3 - k3*p2);
    res(1) = p0*(k3*k4 - k4*p3 + 6*k5) + p1*(-k3*k3 + k3*p3 - k5*k5 + k5*p5 + 12) - p2*(6*k3 - k4*k5 + k4*p5) +
	p3*(k0*k4 + k1*k3 - k1*p3 + 6*k2 - k4*p0) + p4*(-2*k0*k3 + k0*p3 - 2*k2*k5 + k2*p5 + k3*p0 + k5*p2) -
	p5*(6*k0 - k1*k5 + k1*p5 - k2*k4 + k4*p2);
    res(2) = -p0*(-k3*k5 + 6*k4 + k5*p3) + p1*(6*k3 + k4*k5 - k5*p4) + p2*(-k3*k3 + k3*p3 - k4*k4 + k4*p4 + 12) -
	p3*(-k0*k5 + 6*k1 - k2*k3 + k2*p3 + k5*p0) + p4*(6*k0 + k1*k5 + k2*k4 - k2*p4 - k5*p1) +
	p5*(-2*k0*k3 + k0*p3 - 2*k1*k4 + k1*p4 + k3*p0 + k4*p1);
    res(3) = p3*(-k4*k4 + k4*p4 - k5*k5 + k5*p5 + 12) - p4*(-k3*k4 + k3*p4 + 6*k5) + p5*(k3*k5 - k3*p5 + 6*k4);
    res(4) = p3*(k3*k4 - k4*p3 + 6*k5) + p4*(-k3*k3 + k3*p3 - k5*k5 + k5*p5 + 12) - p5*(6*k3 - k4*k5 + k4*p5);
    res(5) = -p3*(-k3*k5 + 6*k4 + k5*p3) + p4*(6*k3 + k4*k5 - k5*p4) + p5*(-k3*k3 + k3*p3 - k4*k4 + k4*p4 + 12);

    res /= 12;
    return;
}

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
