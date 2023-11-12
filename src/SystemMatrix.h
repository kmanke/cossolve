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
#include "MatrixOps.h"
#include "Indexer.h"

#include <eigen3/Eigen/Eigen>

#include <iostream>

namespace cossolve {

enum class BlockDim { row, col };

template <typename MatrixType>
class SystemMatrix
{    
public:
    SystemMatrix(int firstBlockRows, int firstBlockCols);
    ~SystemMatrix();

    // Getter functions
    inline const MatrixType& getMat() const { return mat; }
    inline MatrixType& getMatRef() { return mat; }

    // Returns a reference to the sub matrix at `blockRow`, `blockCol`
    const auto getBlock(int blockRow, int blockCol) const
    { return mat.block(subRows[blockRow].first, subCols[blockCol].first,
		       subRows[blockRow].count, subCols[blockCol].count); }

    // Clears all data in the specified block. Zero entries are not automatically
    // pruned.
    void clearBlock(int blockRow, int blockCol);

    // Clears all the data in the specified dimension
    template <BlockDim dim>
    void clearDim(int blockIndex);

    // Copies all data from `src` to the specified block.
    void copyToBlock(const MatrixType& src, int blockRow, int blockCol);
    
    // Adds a block of size `size` to the matrix. Its index
    // will be equal to the value of `nBlocks<dim>()` before the function call.
    template <BlockDim dim> 
    void addBlock(int size);

    // Adds `count` rows or columns to the specified block.
    // `blockIndex` is the index of the block.
    // `afterIndex` is the index (relative to the start of the block) after which
    // the rows or columns should be inserted.
    template <BlockDim dim>
    void enlargeBlock(int blockIndex, int afterIndex, int count);

    // Returns the size of the `blockIndex`th block along `dim`.
    template <BlockDim dim>
    int blockSize(int blockIndex) const;

    // Adds the values in `src` to the specified block.
    void addToBlock(const MatrixType& src, int blockRow, int blockCol);
    void subtractFromBlock(const MatrixType& src, int blockRow, int blockCol);

    // Returns a reference to a single coefficient within the specified block
    ScalarType& coeffRef(int blockRow, int blockCol, int row, int col);

    // Returns the specified coefficient within the specified block by value.
    ScalarType coeff(int blockRow, int blockCol, int row, int col);

    template <BlockDim dim>
    int getIndex(int blockIndex, int index, int stride = 1, int offset = 0);

    template <unsigned rowStride, unsigned colStride>
    Indexer<rowStride, colStride> getIndexer(int blockRow, int blockCol);
    
private:
    // This structure holds the offsets and sizes of each sub dimension
    struct SubDim
    {
	int first;
	int count;
	int last;
    };
    using SubDimList = std::vector<SubDim>;
    SubDimList subCols;
    SubDimList subRows;

    // Some shared data structures which we make use of
    //const SolverParameters& params;
    //const ConstraintList& fixedConstraints;
    
    // The full system matrix
    MatrixType mat;

    // This function checks for negative indices and adjusts them
    // to point to the correct index.
    inline void adjustIndex(int& blockRow, int& blockCol)
    {
	if (blockRow < 0) { blockRow = subRows.size() + blockRow; }
	if (blockCol < 0) { blockCol = subCols.size() + blockCol; }
	return;
    }
    template <BlockDim dim>
    inline void adjustIndex(int& blockIndex)
    {
	if (blockIndex < 0)
	{
	    if constexpr (dim == BlockDim::row)
	    {
		blockIndex = subRows.size() + blockIndex;
	    }
	    else
	    {
		blockIndex = subCols.size() + blockIndex;
	    }
	}
	return;
    }
};

// Definition of templated functions
template <typename MatrixType>
SystemMatrix<MatrixType>::SystemMatrix(int firstBlockRows, int firstBlockCols)
{
    // Allocate the rows and add the first block to the lists
    mat = MatrixType(firstBlockRows, firstBlockCols);
    SubDim firstRow, firstCol;
    firstRow.first = 0;
    firstRow.count = firstBlockRows;
    firstRow.last = firstBlockRows - 1;
    firstCol.first = 0;
    firstCol.count = firstBlockCols;
    firstCol.last = firstBlockCols - 1;

    subRows.emplace_back(std::move(firstRow));
    subCols.emplace_back(std::move(firstCol));
    
    return;
}

template <typename MatrixType>
SystemMatrix<MatrixType>::~SystemMatrix() { }

template <typename MatrixType>
void SystemMatrix<MatrixType>::clearBlock(int blockRow, int blockCol)
{
    adjustIndex(blockRow, blockCol);
    cossolve::clearBlock(mat, subRows[blockRow].first, subCols[blockCol].first,
			 subRows[blockRow].count, subCols[blockCol].count);
    return;
}

template <typename MatrixType>
template <BlockDim dim>
void SystemMatrix<MatrixType>::clearDim(int blockIndex)
{
    adjustIndex<dim>(blockIndex);
    if constexpr(dim == BlockDim::row)
    {
	for (int i = 0; i < subCols.size(); i++)
	{
	    clearBlock(blockIndex, i);
	}
    }
    else
    {
	for (int i = 0; i < subRows.size(); i++)
	{
	    clearBlock(i, blockIndex);
	}
    }
    return;
}

template <typename MatrixType>
void SystemMatrix<MatrixType>::copyToBlock(const MatrixType& src, int blockRow, int blockCol)
{
    if (blockRow < 0) { blockRow = subRows.size() + blockRow; }
    if (blockCol < 0) { blockCol = subCols.size() + blockCol; }
    cossolve::copyBlock(src, mat, 0, 0, subRows[blockRow].first, subCols[blockCol].first,
			src.rows(), src.cols());
    return;
}

template <typename MatrixType>
void SystemMatrix<MatrixType>::addToBlock(const MatrixType& src, int blockRow, int blockCol)
{
    if (blockRow < 0) { blockRow = subRows.size() + blockRow; }
    if (blockCol < 0) { blockCol = subCols.size() + blockCol; }
    cossolve::addBlock(src, mat, 0, 0, subRows[blockRow].first, subCols[blockCol].first,
		       src.rows(), src.cols());
    return;
}

template <typename MatrixType>
void SystemMatrix<MatrixType>::subtractFromBlock(const MatrixType& src, int blockRow, int blockCol)
{
    if (blockRow < 0) { blockRow = subRows.size() + blockRow; }
    if (blockCol < 0) { blockCol = subCols.size() + blockCol; }
    cossolve::subtractBlock(src, mat, 0, 0, subRows[blockRow].first, subCols[blockCol].first,
			    src.rows(), src.cols());
    return;
}

template <typename MatrixType>
ScalarType& SystemMatrix<MatrixType>::coeffRef(int blockRow, int blockCol, int row, int col)
{
    return mat.coeffRef(getIndex<BlockDim::row>(blockRow, row),
			getIndex<BlockDim::col>(blockCol, col));
}

template <typename MatrixType>
ScalarType SystemMatrix<MatrixType>::coeff(int blockRow, int blockCol, int row, int col)
{
    return mat.coeff(subRows[blockRow].first + row, subCols[blockCol].first + col);
}

template <typename MatrixType>
template <BlockDim dim>
void SystemMatrix<MatrixType>::addBlock(int size)
{
    if constexpr (dim == BlockDim::row)
    {
	SubDim newRow;
	newRow.count = size;
	newRow.first = subRows.back().last + 1;
	newRow.last = newRow.first + newRow.count - 1;
	addDim<AddDimType::row>(mat, newRow.first - 1, newRow.count);
	subRows.emplace_back(std::move(newRow));
    }
    else
    {
	SubDim newCol;
	newCol.count = size;
	newCol.first = subCols.back().last + 1;
	newCol.last = newCol.first + newCol.count - 1;
	addDim<AddDimType::col>(mat, newCol.first - 1, newCol.count);
	subCols.emplace_back(std::move(newCol));
    }
    return;
}

template <typename MatrixType>
template <BlockDim dim>
void SystemMatrix<MatrixType>::enlargeBlock(int blockIndex, int afterIndex, int count)
{
    if constexpr (dim == BlockDim::row)
    {
	if (blockIndex < 0) { blockIndex = subRows.size() + blockIndex; }
	addDim<AddDimType::row>(mat, subRows[blockIndex].first + afterIndex, count);
	subRows[blockIndex].count += count;
	subRows[blockIndex].last += count;
	for (int i = blockIndex + 1; i < subRows.size(); i++)
	{
	    subRows[i].first += count;
	    subRows[i].last += count;
	}
    }
    else
    {
	if (blockIndex < 0) { blockIndex = subCols.size() + blockIndex; }
	addDim<AddDimType::col>(mat, subCols[blockIndex].first + afterIndex, count);
	subCols[blockIndex].count += count;
	subCols[blockIndex].last += count;
	for (int i = blockIndex + 1; i < subCols.size(); i++)
	{
	    subCols[i].first += count;
	    subCols[i].last += count;
	}
    }
    return;
}

template <typename MatrixType>
template <BlockDim dim>
int SystemMatrix<MatrixType>::blockSize(int blockIndex) const
{
    if constexpr(dim == BlockDim::row)
    {
	if (blockIndex < 0) { blockIndex = subRows.size() + blockIndex; }
	return subRows[blockIndex].count;
    }
    else
    {
	if (blockIndex < 0) { blockIndex = subCols.size() + blockIndex; }
	return subCols[blockIndex].count;
    }
}

template <typename MatrixType>
template <BlockDim dim>
int SystemMatrix<MatrixType>::getIndex(int blockIndex, int index, int stride, int offset)
{
    if constexpr(dim == BlockDim::row)
    {
	if (blockIndex < 0) { blockIndex = subRows.size() + blockIndex; }
	return cossolve::getIndex(index, stride, subRows[blockIndex].first + offset,
				  subRows[blockIndex].last + 1);
    }
    else
    {
	if (blockIndex < 0) { blockIndex = subCols.size() + blockIndex; }
	return cossolve::getIndex(index, stride, subCols[blockIndex].first + offset,
				  subRows[blockIndex].last + 1);
    }
}

template <typename MatrixType>
template <unsigned rowStride, unsigned colStride>
Indexer<rowStride, colStride> SystemMatrix<MatrixType>::getIndexer(int blockRow, int blockCol)
{
    return Indexer<rowStride, colStride>(subRows[blockRow].count, subCols[blockCol].count,
					 subRows[blockRow].first, subCols[blockCol].count);
}

} // namespace cossolve

#endif // COSSOLV_SYSTEM_MATRIX_H
