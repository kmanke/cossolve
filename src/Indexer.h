/* This utility class provides helper objects to generate
 * indices with custom stride and offset.
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

#ifndef COSSOLVE_INDEXER_H
#define COSSOLVE_INDEXER_H

#include <eigen3/Eigen/Eigen>

namespace cossolve {

template <unsigned rowStride, unsigned colStride>
class Indexer
{
public:
    Indexer(int nRows, int nCols)
	: nRows(nRows), nCols(nCols), firstRow(0), firstCol(0),
		lastRow(nRows - 1), lastCol(nCols - 1) { }
    Indexer(int nRows, int nCols, int firstRow, int firstCol)
	: nRows(nRows), nCols(nCols), firstRow(firstRow), firstCol(firstCol),
	  lastRow(firstRow + nRows - 1), lastCol(firstCol + nCols - 1) { }								
    ~Indexer() { }

    // Returns the row of the corresponding index.
    inline int row(int index) const
    {
	return firstRow + ((index >= 0) ? (index*rowStride) : (nRows + index*rowStride));
    }
    inline int row(int index, int offset) const
    {
	return row(index) + offset;
    }

    // Returns the column of the corresponding index.
    inline int col(int index) const
    {
	return firstCol + ((index >= 0) ? (index*colStride) : (nCols + index*colStride));
    }
    inline int col(int index, int offset) const
    {
	return col(index) + offset;
    }

    // Updates the size of the indexed matrix. This must be called
    // when the matrix is resized if negative indexing is used.
    inline void resize(int nRows, int nCols)
    {
	this->nRows = nRows;
	this->nCols = nCols;
	return;
    }

    // Updates the base of the indexed matrix.
    inline void rebase(int firstRow, int firstCol)
    {
	this->firstRow = firstRow;
	this->firstCol = firstCol;
	return;
    }

    // Returns a block from the specified matrix
    // The block size is equal to the stride sizes
    template <typename Derived>
    inline Eigen::Block<Derived, rowStride, colStride> block(Eigen::EigenBase<Derived>& mat,
							     int rowIndex, int colIndex)
    {
	return Eigen::Block<Derived, rowStride, colStride>(mat.derived(), row(rowIndex), col(colIndex));
    }
    
    template <typename Derived>
    inline const Eigen::Block<const Derived, rowStride, colStride> block(const Eigen::EigenBase<Derived>& mat,
									 int rowIndex, int colIndex)
    {
	return Eigen::Block<const Derived, rowStride, colStride>(mat.derived(), row(rowIndex), col(colIndex));
    }
    
    // In this overload, the block size is specified
    template <typename Derived>
    inline Eigen::Block<Derived> block(Eigen::EigenBase<Derived>& mat,
				       int rowIndex, int colIndex, int nRows, int nCols)
    {
	return Eigen::Block<Derived>(mat.derived(), row(rowIndex), col(colIndex), nRows, nCols);
    }
    
    template <typename Derived>
    inline const Eigen::Block<const Derived> block(const Eigen::EigenBase<Derived>& mat,
						   int rowIndex, int colIndex, int nRows, int nCols)
    {
	return Eigen::Block<const Derived>(mat.derived(), row(rowIndex), col(colIndex), nRows, nCols);
    }

private:
    int nRows, nCols;
    int firstRow, firstCol;
    int lastRow, lastCol;
};

} // namespace cossolve

#endif // COSSOLVE_INDEX_H
