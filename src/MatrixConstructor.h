/* This class defines operations to construct sparse matrices which
 * feature repeating patterns of coefficients.
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

#ifndef COSSOLVE_MATRIX_CONSTRUCTOR_H
#define COSSOLVE_MATRIX_CONSTRUCTOR_H

#include "config.h"
#include "CossolveTypes.h"

#include <vector>
#include <memory>

namespace cossolve {

class MatrixConstructor
{
public:
    MatrixConstructor() { };
    virtual ~MatrixConstructor() { };

    virtual void construct(SparseType& mat) const = 0;
};

using MatrixConstructorList = std::vector<std::unique_ptr<MatrixConstructor>>;

// Subclass for initializing values along a diagonal
class DiagConstructor : public MatrixConstructor
{
public:
    DiagConstructor(ScalarType value, int startRow, int startCol, int endRow, int stride)
	: value(value), startRow(startRow), startCol(startCol), endRow(endRow), stride(stride) { }
    ~DiagConstructor() { }

    void construct(SparseType& mat) const;

private:
    ScalarType value;
    int startRow;
    int startCol;
    int endRow;
    int stride;
};

// Subclass for initializing values along a row
class RowConstructor : public MatrixConstructor
{
public:
RowConstructor(ScalarType value, int row, int startCol, int endCol, int stride)
    : value(value), row(row), startCol(startCol), endCol(endCol), stride(stride) { }
    ~RowConstructor() { }

    void construct(SparseType& mat) const;

private:
    ScalarType value;
    int row;
    int startCol;
    int endCol;
    int stride;
};

} // namespace cossolve
    
#endif // COSSOLVE_MATRIX_CONSTRUCTOR_H
