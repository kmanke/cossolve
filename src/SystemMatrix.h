/* This class represents the main system matrix for the solver with its
 * submatrices.
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

#ifndef COSSOLVE_SYSTEM_MATRIX_H
#define COSSOLVE_SYSTEM_MATRIX_H

#include "config.h"
#include "CossolveTypes.h"
#include "ConstraintList.h"
#include "SolverParameters.h"

#include <eigen3/Eigen/Eigen>

namespace cossolve {

struct DerivativeMatrixCoefficients
{
    ScalarType firstNodeCurrentCoeff;
    ScalarType firstNodeNextCoeff;
    ScalarType innerNodePrevCoeff;
    ScalarType innerNodeCurrentCoeff;
    ScalarType innerNodeNextCoeff;
    ScalarType lastNodePrevCoeff;
    ScalarType lastNodeCurrentCoeff;
};

struct IntegralMatrixCoefficients
{
    ScalarType firstNodeCurrentCoeff;
    ScalarType firstNodeNextCoeff;
    ScalarType innerNodeFirstCoeff;
    ScalarType innerNodeIntermediateCoeff;
    ScalarType innerNodePrevCoeff;
    ScalarType innerNodeCurrentCoeff;
    ScalarType innerNodeNextCoeff;
    ScalarType lastNodeFirstCoeff;
    ScalarType lastNodeIntermediateCoeff;
    ScalarType lastNodePrevCoeff;
    ScalarType lastNodeCurrentCoeff;
};
    
class SystemMatrix
{    
public:
    SystemMatrix(const SolverParameters& params, const ConstraintList& fixedConstraints);
    ~SystemMatrix();

    // Getter functions
    inline const MatrixType& getMat() const { return mat; }
    inline const MatrixType& getEinv() const { return Einv; }
    inline const MatrixType& getAfK() const { return AfK; }
    
    // Initialization functions
    // Initializes the strain matrix.
    void initStrainMatrix(Eigen::Ref<const TwistType> stiffnessDiag,
			  const DerivativeMatrixCoefficients& derivCoeffs,
			  const IntegralMatrixCoefficients& intCoeffs,
			  SparseLUSolver<MatrixType>& solver);

    // Updates the submatrix dimensions to accomodate an added constraint
    // Note that the matrices are not automatically populated. A call to
    // `updateFixedconstraints` is necessary to do so.
    void addFixedConstraint(int index);

    // Removes the fixed constraint corresponding to the passed
    // constraint number.
    void removeFixedConstraint(int index); // not implemented yet

private:
    // Enum for templated submatrix operations
    enum SubMatrix
    {
	system,
	strain,
	fixedConstraintForce,
	fixedConstraintLocation,
	fixedConstraintActive
    };

    // Some shared data structures which we make use of
    const SolverParameters& params;
    const ConstraintList& fixedConstraints;
    
    // The full system matrix
    MatrixType mat;

    // We need to store copies of some matrices so we can regenerate efficiently
    MatrixType K; // Stiffness
    MatrixType D; // Derivative
    MatrixType E; // Integral
    MatrixType Einv; // E^-1
    MatrixType EinvKD; // E^-1*K*D
    MatrixType Af; // Adjoint strain
    MatrixType AfK; // Af*K

    // Updates the fixed constraint submatrix to reflect changes.
    void updateFixedConstraints();

    // Convenience wrapper around clearBlock which clears a submatrix.
    template <SubMatrix sub, bool prune>
    void clearBlock();
  
    // Helper functions for indexing the various matrix types
    template <SubMatrix sub>
    int firstRow();
    template <SubMatrix sub>
    int firstCol();
    template <SubMatrix sub>
    int lastRow();
    template <SubMatrix sub>
    int lastCol();
    template <SubMatrix sub>
    int nRows();
    template <SubMatrix sub>
    int nCols();
    template <SubMatrix sub>
    int rowIndex(int row, int stride = 1);
    template <SubMatrix sub>
    int colIndex(int col, int stride = 1);
    
    // Internal functions
    void initStiffnessMatrix(Eigen::Ref<const TwistType> diag);
    void initDerivativeMatrix(const DerivativeMatrixCoefficients& coeffs);
    void initIntegralMatrix(const IntegralMatrixCoefficients& coeffs);
};

// Definition of templated functions
template <SystemMatrix::SubMatrix sub>
int SystemMatrix::firstRow()
{
    if constexpr (sub == system)
    {
	return 0;
    }
    if constexpr (sub == strain)
    {
	return 0;
    }
    else if constexpr (sub == fixedConstraintForce)
    {
	return 0;
    }
    else if constexpr (sub == fixedConstraintLocation)
    {
	return lastRow<strain>() + 1;
    }
    else if constexpr (sub == fixedConstraintActive)
    {
	return lastRow<fixedConstraintForce>() + 1;
    }
    else
    {
	return 0;
    }
}
    
template <SystemMatrix::SubMatrix sub>
int SystemMatrix::lastRow()
{
    return firstRow<sub>() + nRows<sub>() - 1;
}

template <SystemMatrix::SubMatrix sub>
int SystemMatrix::firstCol()
{
    if constexpr (sub == system)
    {
	return 0;
    }
    if constexpr (sub == strain)
    {
	return 0;
    }
    else if constexpr (sub == fixedConstraintForce)
    {
	return lastCol<strain>() + 1;
    }
    else if constexpr (sub == fixedConstraintLocation)
    {
	return 0;
    }
    else if constexpr (sub == fixedConstraintActive)
    {
	return lastCol<fixedConstraintLocation>() + 1;
    }
    else
    {
	return 0;
    }
}

template <SystemMatrix::SubMatrix sub>
int SystemMatrix::lastCol()
{
    return firstCol<sub>() + nCols<sub>() - 1;
}
    
template <SystemMatrix::SubMatrix sub>
int SystemMatrix::nRows()
{
    if constexpr (sub == system)
    {
	return nRows<strain>() + nRows<fixedConstraintLocation>();
    }
    if constexpr (sub == strain)
    {
	return params.nNodes() * twistLength;
    }
    else if constexpr (sub == fixedConstraintForce)
    {
	return params.nNodes() * twistLength;
    }
    else if constexpr (sub == fixedConstraintLocation)
    {
	return fixedConstraints.size();
    }
    else if constexpr (sub == fixedConstraintActive)
    {
	return fixedConstraints.size();
    }
    else
    {
	return 0;
    }
}
    
template <SystemMatrix::SubMatrix sub>
int SystemMatrix::nCols()
{
    if constexpr (sub == system)
    {
	return nCols<strain>() + nCols<fixedConstraintForce>();
    }
    else if constexpr (sub == strain)
    {
	return params.nNodes() * twistLength;
    }
    else if constexpr (sub == fixedConstraintForce)
    {
	return fixedConstraints.size();
    }
    else if constexpr (sub == fixedConstraintLocation)
    {
	return params.nNodes() * twistLength;
    }
    else if constexpr (sub == fixedConstraintActive)
    {
	return fixedConstraints.size();
    }
    else
    {
	return 0;
    }
}

template <SystemMatrix::SubMatrix sub>
int SystemMatrix::rowIndex(int row, int stride)
{
    return (row >= 0) ? (firstRow<sub>() + row*stride) : (lastRow<sub>() + 1 + row*stride);
}

template <SystemMatrix::SubMatrix sub>
int SystemMatrix::colIndex(int col, int stride)
{
    return (col >= 0) ? (firstCol<sub>() + col*stride) : (lastCol<sub>() + 1 + col*stride);
}

} // namespace cossolve

#endif // COSSOLV_SYSTEM_MATRIX_H
