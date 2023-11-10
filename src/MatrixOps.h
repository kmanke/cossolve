/* This file and its associated .cpp file implement some useful
 * matrix operations which are not natively included in Eigen.
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

#ifndef COSSOLVE_MATRIX_OPS_H
#define COSSOLVE_MATRIX_OPS_H

#include "CossolveTypes.h"
#include "config.h"

#include <vector>

namespace cossolve {

// Returns a block reference into the specified vector.
Eigen::Ref<VectorType> vectorRef(Eigen::Ref<VectorType> vec, int index, int stride = 1);
Eigen::Ref<const VectorType> vectorRefConst(Eigen::Ref<const VectorType> vec, int index, int stride = 1);

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

// Adds `nDims` starting after `afterDim` (which can be -1 to insert at the top / left)
enum AddDimType { addRow, addCol };
template <AddDimType dim>
void addDim(MatrixType& src, int afterDim, int nDims);

// Implementation of templated functions follows.
template <bool prune>
void clearBlock(MatrixType& dst, int startRow, int startCol, int nRows, int nCols)
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
void copyBlock(const MatrixType& src, MatrixType& dst,
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

template <AddDimType dim>
void addDim(MatrixType& dst, int afterDim, int nDims)
{
    // Start by resizing the matrix to add the rows
    if constexpr (dim == addRow)
    {
	dst.conservativeResize(dst.rows() + nDims, dst.cols());
    }
    else
    {
	dst.conservativeResize(dst.rows(), dst.cols() + nDims);
    }
    // If the we're adding to the inner dimension, we simply add nRows to any inner indices which are
    // greater than `afterDim`.
    if constexpr ((!MatrixType::IsRowMajor && dim == addRow) ||
		  (MatrixType::IsRowMajor && dim == addCol))
    {
	MatrixType::StorageIndex* innerPtr = dst.innerIndexPtr();
	// The one-past-the-end outer start points to one-past-the-end of the inner indices.
	MatrixType::StorageIndex* innerEnd = innerPtr + dst.outerIndexPtr()[dst.outerSize()];
	while (innerPtr < innerEnd)
	{
	    if (*innerPtr > afterDim)
	    {
		*innerPtr += nDims;
	    }
	    innerPtr++;
	}
    }
    // Otherwise, we're adding to the outer, and need to repoint each outer, starting from
    // the back.
    else
    {
	MatrixType::StorageIndex* outerPtr = dst.outerIndexPtr();
	MatrixType::StorageIndex* nnzPtr = dst.innerNonZeroPtr();
	// First, repoint the outers after the inserted dims
	for (int i = dst.outerSize() - 1; i > afterDim + nDims; i--)
	{
	    outerPtr[i] = outerPtr[i - nDims];
	    if (nnzPtr)
	    {
		nnzPtr[i] = nnzPtr[i - nDims];
	    }
	}
	// Now, the inserted outers all just point to the same spot (until something gets inserted)
	for (int i = afterDim + 1; i <= afterDim + nDims; i++)
	{
	    outerPtr[i] = outerPtr[afterDim + nDims + 1];
	    if (nnzPtr)
	    {
		    nnzPtr[i] = 0;
	    }
	}
    }
    return;
}

} // namespace cossolve
    
#endif // COSSOLVE_MATRIX_OPS_H
