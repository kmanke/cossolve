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

class MatrixConstructorList
{
public:
    MatrixConstructorList() { }
    ~MatrixConstructorList() { }

    // Adds a new MatrixConstructor to the list by calling the constructor
    // with the specified args.
    template <typename MatrixConstructorType, typename ...Args>
    void add(Args... args);

    // Clears the list.
    inline void clear() { list.clear(); }

    // Iterator access
    inline auto begin() { return list.begin(); }
    inline auto end() { return list.end(); }
    inline auto cbegin() const { return list.cbegin(); }
    inline auto cend() const { return list.cend(); }
    
private:
    std::vector<std::unique_ptr<MatrixConstructor>> list;
};

// Implementation of templated member functions
template <typename MatrixConstructorType, typename ...Args>
void MatrixConstructorList::add(Args... args)
{
    list.emplace_back(
	std::make_unique<MatrixConstructorType>(std::forward<Args>(args)...));
    return;
}

} // namespace cossolve
    
#endif // COSSOLVE_MATRIX_CONSTRUCTOR_H
