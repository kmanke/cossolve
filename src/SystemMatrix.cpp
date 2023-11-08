/* Implementation of the SystemMatrix class.
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

#include "SystemMatrix.h"

namespace cossolve {

SystemMatrix::SystemMatrix()
{
}

SystemMatrix::~SystemMatrix()
{
}

template <bool prune>
void SystemMatrix::clearBlock(MatrixType& dst, int startRow, int startCol, int nRows, int nCols)
{
    // Iterate through and zero all values
    int endRow = startRow + nRows;
    int endCol = startCol + nCols;
    for (int outer = 0; outer < dst.outerSize(); outer++)
    {
	auto innerIt = MatrixType::InnerIterator(dst, outer);
	// Check if the outer is in our range
	if ((MatrixType::IsRowMajor && innerIt.row() >= startRow && innerIt.row() < endRow) ||
	    (!MatrixType::IsRowMajor && innerIt.col() >= startCol && innerIt.col() < endCol))
	{
	    for (; innerIt; ++innerIt)
	    {
		// Check if we're past the first row / col
		if ((MatrixType::IsRowMajor && innerIt.col() >= startCol) ||
		    (!MatrixType::IsRowMajor && innerIt.row() >= startRow))
		{
		    // Check if we're past the last row / col
		    // Since inners are guaranteed to be ascending, we can break here
		    if ((MatrixType::IsRowMajor && innerIt.col() >= endCol) ||
			(!MatrixType::IsRowMajor && innerIt.row() >= endRow))
		    {
			break;
		    }
		    // Found a value in the block, zero it out
		    innerIt.valueRef() = 0;
		}
	    }
	}
    }
    // Prune if that version of the function is called
    if constexpr (prune)
    {
	auto keepFunc = [&](const Eigen::Index& row, const Eigen::Index& col, const MatrixType::Scalar& value)
	    {
		if ((row >= startRow && row < endRow) &&
		    (col >= startCol && col < endCol) && value == 0)
		{
		    return true;
		}
		return false;
	    };
	dst.prune(keepFunc);
    }
    return;
}

template <bool prune>
void SystemMatrix::copyBlock(const MatrixType& src, MatrixType& dst,
			     int srcStartRow, int srcStartCol, int dstStartRow, int dstStartCol,
			     int nRows, int nCols)
{
    // Iterate through and copy all values
    int endRow = srcStartRow + nRows;
    int endCol = srcStartCol + nCols;
    // First, clear the destination block
    clearBlock<false>(dst, dstStartRow, dstStartCol, nRows, nCols);
    for (int outer = 0; outer < src.outerSize(); outer++)
    {
	auto innerIt = MatrixType::InnerIterator(src, outer);
	// Check if the outer is in our range
	if ((MatrixType::IsRowMajor && innerIt.row() >= srcStartRow && innerIt.row() < endRow) ||
	    (!MatrixType::IsRowMajor && innerIt.col() >= srcStartCol && innerIt.col() < endCol))
	{
	    for (; innerIt; ++innerIt)
	    {
		// Check if we're past the first row / col
		if ((MatrixType::IsRowMajor && innerIt.col() >= srcStartCol) ||
		    (!MatrixType::IsRowMajor && innerIt.row() >= srcStartRow))
		{
		    // Check if we're past the last row / col
		    // Since inners are guaranteed to be ascending, we can break here
		    if ((MatrixType::IsRowMajor && innerIt.col() >= endCol) ||
			(!MatrixType::IsRowMajor && innerIt.row() >= endRow))
		    {
			break;
		    }
		    // Found a value in the src block. Copy it to dst.
		    dst.coeffRef(innerIt.row(), innerIt.col()) = innerIt.value();
		}
	    }
	}
    }
    // Prune if that version of the function is called
    if constexpr (prune)
    {
	int dstEndRow = dstStartRow + nRows;
	int dstEndCol = dstStartCol + nCols;
	auto keepFunc = [&](const Eigen::Index& row, const Eigen::Index& col, const MatrixType::Scalar& value)
	    {
		if ((row >= dstStartRow && row < dstEndRow) &&
		    (col >= dstStartCol && col < dstEndCol))
		{
		    return false;
		}
		return true;
	    };
	dst.prune(keepFunc);
    }
    return;    
}   

void SystemMatrix::initStiffnessMatrix(Eigen::Ref<const TwistType> diag)
{
    // Reinitializing the stiffness matrix clears the strain matrix
    clearBlock<true>(mat, firstRow<strain>(), firstCol<strain>(), nRows<strain>(), nCols<strain>());
    for (int i = firstCol<strain>(); i <= lastCol<strain>(); i++)
    {
	mat.coeffRef(i, i) = diag(i % twistLength);
    }
    
    return;
}

void SystemMatrix::initDerivativeMatrix(const DerivativeMatrixCoefficients& coeffs)
{
    // We make a temporary copy before multiplying it with the stiffness matrix
    MatrixType derivativeMatrix(nRows<strain>(), nCols<strain>());
    int i = firstRow<strain>();
    // First node
    for (; i < rowIndex<strain>(1, twistLength); i++)
    {
	mat.coeffRef(i, i) = coeffs.firstNodeCurrentCoeff;
	mat.coeffRef(i, i + twistLength) = coeffs.firstNodeNextCoeff;
    }
    // Inner nodes
    for (; i < rowIndex<strain>(-1, twistLength); i++)
    {
	mat.coeffRef(i, i - twistLength) = coeffs.innerNodePrevCoeff;
	mat.coeffRef(i, i) = coeffs.innerNodeCurrentCoeff;
	mat.coeffRef(i, i + twistLength) = coeffs.innerNodeNextCoeff;
    }
    // Last node
    for (; i <= lastRow<strain>(); i++)
    {
	mat.coeffRef(i, i - twistLength) = coeffs.lastNodePrevCoeff;
	mat.coeffRef(i, i) = coeffs.lastNodeCurrentCoeff;
    }
    // Now we multiply the stiffness and derivative matrices together
    derivativeMatrix = mat.block(firstRow<strain>(), firstCol<strain>(),
				 nRows<strain>(), nCols<strain>()) * derivativeMatrix;
    copyBlock<true>(derivativeMatrix, mat, 0, 0, firstRow<strain>(), firstCol<strain>(),
		    nRows<strain>(), nCols<strain>());
    return;
}

void SystemMatrix::initIntegralMatrix(const IntegralMatrixCoefficients& coeffs)
{
    // We make a temporary copy before multiplying it with the stiffness matrix
    MatrixType integralMatrix(nRows<strain>(), nCols<strain>());
    int i = firstRow<strain>();
    // First node
    for (; i < rowIndex<strain>(1, twistLength); i++)
    {
	mat.coeffRef(i, i) = coeffs.firstNodeCurrentCoeff;
	mat.coeffRef(i, i + twistLength) = coeffs.firstNodeNextCoeff;
    }
    // Inner nodes
    for (; i < rowIndex<strain>(-1, twistLength); i++)
    {
	int j = firstCol<strain>() + i % twistLength;
	mat.coeffRef(i, j) = coeffs.innerNodeFirstCoeff;
	j += twistLength;
	for (; j < i - twistLength; j += twistLength)
	{
	    mat.coeffRef(i, j) = coeffs.innerNodeIntermediateCoeff;
	}
	mat.coeffRef(i, i - twistLength) = coeffs.innerNodePrevCoeff;
	mat.coeffRef(i, i) = coeffs.innerNodeCurrentCoeff;
	mat.coeffRef(i, i + twistLength) = coeffs.innerNodeNextCoeff;
    }
    // Last node
    for (; i <= lastRow<strain>(); i++)
    {
	int j = firstCol<strain>() + i % twistLength;
	mat.coeffRef(i, j) = coeffs.lastNodeFirstCoeff;
	j += twistLength;
	for (; j < i - twistLength; j += twistLength)
	{
	    mat.coeffRef(i, j) = coeffs.lastNodeIntermediateCoeff;
	}
	mat.coeffRef(i, i - twistLength) = coeffs.lastNodePrevCoeff;
	mat.coeffRef(i, i) = coeffs.lastNodeCurrentCoeff;
    }
    // The integral matrix now gets inverted, multiplied with the strain and derivative
    // matrices, and then copied into the system matrix.
    
    return;
}

} // namespace cossolve
