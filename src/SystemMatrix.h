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
//    inline const MatrixType& getEinv() const { return Einv; }
//    inline const MatrixType& getAfK() const { return AfK; }

    // Returns a reference to the sub matrix at `blockRow`, `blockCol`
    const auto getBlock(int blockRow, int blockCol) const
    { return mat.block(subRows[blockRow].first, subCols[blockCol].first,
		       subRows[blockRow].count, subCols[blockCol].count); }

    // Clears all data in the specified block. Zero entries are not automatically
    // pruned.
    void clearBlock(int blockRow, int blockCol);

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

    // Returns a reference to a single coefficient within the specified block
    ScalarType& coeffRef(int blockRow, int blockCol, int rowIndex, int colIndex,
			 int stride = 1, int offset = 0);

    // Returns the specified coefficient within the specified block by value.
    ScalarType coeff(int blockRow, int blockCol, int row, int col);

    template <BlockDim dim>
    int getIndex(int blockIndex, int index, int stride = 1, int offset = 0);
    
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

    // Updates the fixed constraint submatrix to reflect changes.
    void updateFixedConstraints();
    
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

    // We need to store copies of some matrices so we can regenerate efficiently
/*    MatrixType K; // Stiffness
    MatrixType D; // Derivative
    MatrixType E; // Integral
    MatrixType Einv; // E^-1
    MatrixType EinvKD; // E^-1*K*D
    MatrixType Af; // Adjoint strain
    MatrixType AfK; // Af*K */
  
    /*// Helper functions for indexing the various matrix types
    template <SubMatrix sub>
    int firstIndex();
    template <SubMatrix sub>
    int lastIndex();
    template <SubMatrix sub>
    int nDims();
    template <SubMatrix sub>
    int subIndex(int index, int stride = 1); */
    
    // Internal functions
    void initStiffnessMatrix(Eigen::Ref<const TwistType> diag);
    void initDerivativeMatrix(const DerivativeMatrixCoefficients& coeffs);
    void initIntegralMatrix(const IntegralMatrixCoefficients& coeffs);
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
    if (blockRow < 0) { blockRow = subRows.size() + blockRow; }
    if (blockCol < 0) { blockCol = subCols.size() + blockCol; }
    cossolve::clearBlock(mat, subRows[blockRow].first, subCols[blockCol].first,
			 subRows[blockRow].count, subCols[blockCol].count);
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
ScalarType& SystemMatrix<MatrixType>::coeffRef(int blockRow, int blockCol, int row, int col, int stride, int offset)
{
    return mat.coeffRef(getIndex<BlockDim::row>(blockRow, row, stride, offset),
			getIndex<BlockDim::col>(blockCol, col, stride, offset));
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
	newRow.first = subCols.back().last + 1;
	newRow.last = newRow.first + newRow.count - 1;
	subRows.emplace_back(std::move(newRow));
	enlargeBlock<BlockDim::row>(subRows.size() - 1, -1, size);
    }
    else
    {
	SubDim newCol;
	newCol.count = size;
	newCol.first = subCols.back().last + 1;
	newCol.last = newCol.first + newCol.count - 1;
	subCols.emplace_back(std::move(newCol));
	enlargeBlock<BlockDim::col>(subCols.size() - 1, -1, size);
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


/*template <SystemMatrix::SubMatrix sub>
int SystemMatrix::firstIndex()
{
    if constexpr (sub == SubMatrix::system ||
		  sub == SubMatrix::twist)
    {
	return 0;
    }
    else if constexpr (sub == SubMatrix::strain)
    {
	return lastIndex<SubMatrix::twist>() + 1;
    }
    else if constexpr (sub == SubMatrix::fixed)
    {
	return lastIndex<SubMatrix::strain>() + 1;
    }
    else
    {
	return 0;
    }
}
    
template <SystemMatrix::SubMatrix sub>
int SystemMatrix::lastIndex()
{
    return firstIndex<sub>() + nDims<sub>() - 1;
}
    
template <SystemMatrix::SubMatrix sub>
int SystemMatrix::nDims()
{
    if constexpr (sub == SubMatrix::system)
    {
	return nDims<SubMatrix::twist>() + nDims<SubMatrix::strain>() + nDims<SubMatrix::fixed>();
    }
    else if constexpr (sub == SubMatrix::twist || sub == SubMatrix::strain)
    {
	return params.nNodes() * twistLength;
    }
    else if constexpr (sub == SubMatrix::twist)
    {
	return fixedConstraints.size() * twistLength;
    }
    else
    {
	return 0;
    }
}
    
template <SystemMatrix::SubMatrix sub>
int SystemMatrix::subIndex(int index, int stride)
{
    return (index >= 0) ? (firstIndex<sub>() + index*stride) : (lastIndex<sub>() + 1 + index*stride);
    }*/

} // namespace cossolve

#endif // COSSOLV_SYSTEM_MATRIX_H
