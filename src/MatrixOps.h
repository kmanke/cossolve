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

// Convenience function to calculate indices within matrices or vectors.
// The returned value is `index` * `stride` + `offset`.
// If index is less than 0, then the function returns the index relative to
// the end. `endIndex` must be provided in this case.
inline int getIndex(int index, int stride = 1, int offset = 0, int endIndex = 0)
{
    return (index < 0) ? (endIndex + index*stride + offset) : index*stride + offset;
}

// Returns a block reference into the specified vector.
Eigen::Ref<VectorType> vectorRef(Eigen::Ref<VectorType> vec, int index, int stride = 1);
Eigen::Ref<const VectorType> vectorRefConst(Eigen::Ref<const VectorType> vec, int index, int stride = 1);

// Sets all values in the specified block to zero.
template <typename Matrix>
void clearBlock(Matrix& dst, int startRow, int startCol, int nRows, int nCols);

// Prunes all zero-valued entries in the specified block.
void pruneBlock(SparseType& dst, int startRow, int startCol, int nRows, int nCols);

// Sets all values in the specified block to the values in the
// reference matrix.
// If `prune` is true, the block will be pruned after the operation.
template <typename Matrix>
void copyBlock(const Matrix& src, Matrix& dst,
	       int srcStartRow, int srcStartCol, int dstStartRow, int dstStartCol,
	       int nRows, int nCols);

// Adds `nDims` starting after `afterDim` (which can be -1 to insert at the top / left)
enum class AddDimType { row, col };
template <AddDimType dim, typename Matrix>
void addDim(Matrix& src, int afterDim, int nDims);

template <typename Matrix>
void addBlock(const Matrix& src, Matrix& dst,
	      int srcStartRow, int srcStartCol, int dstStartRow, int dstStartCol,
	      int nRows, int nCols);

// Implementation of templated functions follows.

// Anonymous namespace for specializations of addDim
namespace {

template <AddDimType dim>
void addDimSparse(SparseType& dst, int afterDim, int nDims)
{
        // Start by resizing the matrix to add the rows
    if constexpr (dim == AddDimType::row)
    {
	dst.conservativeResize(dst.rows() + nDims, dst.cols());
    }
    else
    {
	dst.conservativeResize(dst.rows(), dst.cols() + nDims);
    }
    // If the we're adding to the inner dimension, we simply add nRows to any inner indices which are
    // greater than `afterDim`.
    if constexpr ((!SparseType::IsRowMajor && dim == AddDimType::row) ||
		  (SparseType::IsRowMajor && dim == AddDimType::col))
    {
	SparseType::StorageIndex* innerPtr = dst.innerIndexPtr();
	// The one-past-the-end outer start points to one-past-the-end of the inner indices.
	SparseType::StorageIndex* innerEnd = innerPtr + dst.outerIndexPtr()[dst.outerSize()];
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
	SparseType::StorageIndex* outerPtr = dst.outerIndexPtr();
	SparseType::StorageIndex* nnzPtr = dst.innerNonZeroPtr();
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

template <AddDimType dim>
void addDimDense(DenseType& dst, int afterDim, int nDims)
{
    if constexpr (dim == AddDimType::row)
    {
	DenseType temp = dst.block(afterDim + 1, 0, dst.rows() - (afterDim + 1), dst.cols());
	dst.conservativeResize(dst.rows() + nDims, dst.cols());
	dst.bottomRows(dst.rows() - (afterDim + 1)) = temp;
    }
    else if constexpr (dim == AddDimType::col)
    {
	DenseType temp = dst.block(0, afterDim + 1, dst.rows(), dst.cols() - (afterDim + 1));
	dst.conservativeResize(dst.rows(), dst.cols() + nDims);
	dst.rightCols(dst.cols() - (afterDim + 1)) = temp;
    }
    return;
}

} // End anonymous namespace

template <AddDimType dim, typename Matrix>
void addDim(Matrix& dst, int afterDim, int nDims)
{
    if constexpr (std::is_same_v<Matrix, SparseType>)
    {
	addDimSparse<dim>(dst, afterDim, nDims);
    }
    else if constexpr (std::is_same_v<Matrix, DenseType>)
    {
	addDimDense<dim>(dst, afterDim, nDims);
    }
    return;
}


} // namespace cossolve
    
#endif // COSSOLVE_MATRIX_OPS_H
