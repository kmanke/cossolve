/* Implementation of MatrixOps.
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

#include "MatrixOps.h"

namespace cossolve {

void constructMatrix(const MatrixConstructorList& constructors, SparseType& mat)
{
    for (auto it = constructors.cbegin(); it != constructors.cend(); ++it)
    {
	(*it)->construct(mat);
    }
    return;
}

Eigen::Ref<VectorType> vectorRef(Eigen::Ref<VectorType> vec, int index, int stride)
{
    return vec.block((index >= 0) ? (index * stride) : (vec.rows() + index * stride), 0,
		     stride, 1);
}

// Specialization of `clearBlock` for sparse matrices.
template <>
void clearBlock<SparseType>(SparseType& dst, int startRow, int startCol, int nRows, int nCols)
{
    // Iterate through and zero all values
    int startOuter = (SparseType::IsRowMajor) ? startRow : startCol;
    int endOuter = startOuter + ((SparseType::IsRowMajor) ? nRows : nCols);
    int startInner = (SparseType::IsRowMajor) ? startCol : startRow;
    int endInner = startInner + ((SparseType::IsRowMajor) ? nCols : nRows);

    for (int outer = startOuter; outer < endOuter; outer++)
    {
	for (auto innerIt = SparseType::InnerIterator(dst, outer); innerIt; ++innerIt)
	{
	    // Check if we're past the first row / col
	    if (innerIt.index() >= startInner)
	    {
		    // Check if we're past the last row / col
		    // Since inners are guaranteed to be ascending, we can break here
		if (innerIt.index() >= endInner)
		{
		    break;
		}
		// Found a value in the block, zero it out
		innerIt.valueRef() = 0;
	    }
	}
    }
    return;
}

// Specialization of `clearBlock` for dense matrices.
template <>
void clearBlock<DenseType>(DenseType& dst, int startRow, int startCol, int nRows, int nCols)
{
    dst.block(startRow, startCol, nRows, nCols).setZero();
    return;
}

// Specialization of `copyBlock` for sparse matrices
template<>
void copyBlock<SparseType>(const SparseType& src, SparseType& dst,
	      int srcStartRow, int srcStartCol, int dstStartRow, int dstStartCol,
	      int nRows, int nCols)
{
    if (srcStartRow < 0) { srcStartRow = src.rows() + srcStartRow; }
    if (srcStartCol < 0) { srcStartCol = src.cols() + srcStartCol; }
    if (dstStartRow < 0) { dstStartRow = dst.rows() + dstStartRow; }
    if (dstStartCol < 0) { dstStartCol = dst.cols() + dstStartCol; }
    int startOuter = (SparseType::IsRowMajor) ? srcStartRow : srcStartCol;
    int endOuter = startOuter + ((SparseType::IsRowMajor) ? nRows : nCols);
    int startInner = (SparseType::IsRowMajor) ? srcStartCol : srcStartRow;
    int endInner = startInner + ((SparseType::IsRowMajor) ? nCols : nRows);
    int deltaRow = dstStartRow - srcStartRow;
    int deltaCol = dstStartCol - srcStartCol;
    clearBlock(dst, dstStartRow, dstStartCol, nRows, nCols);
    for (int outer = startOuter; outer < endOuter; outer++)
    {
	for (auto innerIt = SparseType::InnerIterator(src, outer); innerIt; ++innerIt)
	{
	    // Check if we're past the first row / col
	    if (innerIt.index() >= startInner)
	    {	
		// Check if we're past the last row / col
		// Since inners are guaranteed to be ascending, we can break here
		if (innerIt.index() >= endInner)
		{
		    break;
		}
		// Found a value in the src block. Add it to dst.
		dst.coeffRef(innerIt.row() + deltaRow, innerIt.col() + deltaCol) = innerIt.value();
	    }
	}
    }
    return;
}

// Specialization of `copyblock` for dense matrices.
template<>
void copyBlock<DenseType>(const DenseType& src, DenseType& dst,
	      int srcStartRow, int srcStartCol, int dstStartRow, int dstStartCol,
	      int nRows, int nCols)
{
    if (srcStartRow < 0) { srcStartRow = src.rows() + srcStartRow; }
    if (srcStartCol < 0) { srcStartCol = src.cols() + srcStartCol; }
    if (dstStartRow < 0) { dstStartRow = dst.rows() + dstStartRow; }
    if (dstStartCol < 0) { dstStartCol = dst.cols() + dstStartCol; }
    dst.block(dstStartRow, dstStartCol, nRows, nCols) =
	src.block(srcStartRow, srcStartCol, nRows, nCols);
    return;
}

} // namespace cossolve
