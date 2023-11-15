/* Implementation of MatrixConstructor classes.
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

#include "MatrixConstructor.h"

namespace cossolve {

void DiagConstructor::construct(SparseType& mat) const
{
    int row = (startRow >= 0) ? startRow : (mat.rows() + startRow);
    int col = (startCol >= 0) ? startCol : (mat.cols() + startCol);
    int endTemp = (endRow >= 0) ? endRow : (mat.rows() + endRow);
    while (row < endTemp)
    {
	if (col < mat.cols())
	{
	    mat.coeffRef(row, col) = value;
	    col += stride;
	    row += stride;
	}
	else
	{
	    break;
	}
    }
    return;
}

void RowConstructor::construct(SparseType& mat) const
{
    int col = (startCol >= 0) ? startCol : (mat.cols() + startCol);
    int rowTemp = (row >= 0) ? row : (mat.rows() + row);
    int endTemp = (endCol >= 0) ? endCol  : (mat.cols() + endCol);
    while (col < endTemp)
    {
	mat.coeffRef(rowTemp, col) = value;
	col += stride;
    }
    return;
}

} // namespace cossolve
