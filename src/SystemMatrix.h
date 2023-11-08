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

#include "CossolveTypes.h"

#include <eigen3/Eigen/Eigen>
#include <cassert>

namespace cossolve {

// Convenience overload for adding two indices
std::pair<int, int> operator+(const std::pair<int, int> a, const std::pair<int, int> b)
{
    return std::pair<int, int>(a.first + b.first, a.second + b.second);
}

struct DerivativeMatrixCoefficients
{
    int firstNodeCurrentCoeff;
    int firstNodeNextCoeff;
    int innerNodePrevCoeff;
    int innerNodeCurrentCoeff;
    int innerNodeNextCoeff;
    int lastNodePrevCoeff;
    int lastNodeIntermediateCoeff;
    int lastNodeCurrentCoeff;
};

struct IntegralMatrixCoefficients
{
    int firstNodeCurrentCoeff;
    int firstNodeNextCoeff;
    int innerNodeFirstCoeff;
    int innerNodeIntermediateCoeff;
    int innerNodePrevCoeff;
    int innerNodeCurrentCoeff;
    int innerNodeNextCoeff;
    int lastNodeFirstCoeff;
    int lastNodePrevCoeff;
    int lastNodeIntermediateCoeff;
    int lastNodeCurrentCoeff;
};
    
class SystemMatrix
{
public:
    SystemMatrix();
    ~SystemMatrix();

    // Getter functions
    inline Eigen::Ref<MatrixType> getMat() { return mat; }

    // Initialization functions
    // Initializes the strain matrix.
    void initStrainMatrix(int nNodes, SingleMatrixType nodeStiffness,
			  const PairList& derivativeFirstRow,
			  const PairList& derivativeInnerRow,
			  const PairList& derivativeLastRow,
			  const PairList& integralFirstRow,
			  const PairList& integralInnerRow,
			  const PairList& integralLastRow);

    // Adds a fixed constraint to `node`. `node`'s coordinate transformation
    // will be forced to equal `g`.
    void addFixedConstraint(int node, Eigen::Ref<const CoordType> g);

    // Removes the fixed constraint corresponding to the passed
    // constraint number.
    void removeFixedConstraint(int constraint);

    // Deactivates the fixed constraint. The constraint will still exist in
    // the system matrix, but it will not be enforced.
    void deactivateFixedConstraint(int constraint);

    // Activates the fixed constraint.
    void activateFixedConstraint(int constraint);
    
private:
    // Type aliases
    using Index = std::pair<int, int>;
    enum SubMatrix
    {
	strain,
	fixedConstraintForce,
	fixedConstraintLocation
    };
    
    // The full system matrix
    MatrixType mat;

    // Sizes
    int nNodes;
    int nFixedConstraints;
    int nContactNodes;

    // We need to store copies of some matrices so we can regenerate efficiently
    MatrixType EinvKD;
    MatrixType Af;
    
    // Sets all values in the specified block to zero.
    // If prune is true, the block will be pruned after the operation.
    template <bool prune>
    void clearBlock(MatrixType& dst, int startRow, int startCol, int nRows, int nCols);

    // Sets all values in the specified block to the values in the
    // reference matrix.
    // If `prune` is true, the block will be pruned after the operation.
    template <bool prune>
    void copyBlock(const MatrixType& src, MatrixType& dst,
		   int srcStartRow, int srcStartCol, int dstStartRow, int dstStartCol,
		   int nRows, int nCols);
    
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
    if constexpr (sub == strain)
    {
	return nNodes * twistLength;
    }
    else if constexpr (sub == fixedConstraintForce)
    {
	return nNodes * twistLength;
    }
    else if constexpr (sub == fixedConstraintLocation)
    {
	return nFixedConstraints;
    }
    else
    {
	return 0;
    }
}
    
template <SystemMatrix::SubMatrix sub>
int SystemMatrix::nCols()
{
    if constexpr (sub == strain)
    {
	return nNodes * twistLength;
    }
    if constexpr (sub == fixedConstraintForce)
    {
	return nFixedConstraints;
    }
    if constexpr (sub == fixedConstraintLocation)
    {
	return nNodes * twistLength;
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
