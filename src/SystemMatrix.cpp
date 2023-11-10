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
#include "MatrixOps.h"

#include <iostream>

namespace cossolve {

SystemMatrix::SystemMatrix(const SolverParameters& params,
			   const ConstraintList& fixedConstraints)
    : params(params), fixedConstraints(fixedConstraints)
{
    // Allocate the full system matrix but don't initialize yet
    mat = MatrixType(nDims<SubMatrix::system>(), nDims<SubMatrix::system>());
    return;
}

SystemMatrix::~SystemMatrix() { }

void SystemMatrix::initStrainMatrix(Eigen::Ref<const TwistType> stiffnessDiag,
				    const DerivativeMatrixCoefficients& derivCoeffs,
				    const IntegralMatrixCoefficients& intCoeffs,
				    SparseLUSolver<MatrixType>& solver)
{
    // Start by constructing our basic building blocks
    clearBlock<SubMatrix::twist, SubMatrix::twist, true>();
    initStiffnessMatrix(stiffnessDiag);
    initDerivativeMatrix(derivCoeffs);
    initIntegralMatrix(intCoeffs);

    // The integral matrix now gets inverted, multiplied with the strain and derivative
    // matrices, and then copied into the system matrix.
    MatrixType I = MatrixType(E.rows(), E.cols());
    I.setIdentity();
    solver.compute(E);
    Einv = E;
    Einv = solver.solve(I);
    Einv = Einv.pruned(); // Remove zero elements
    EinvKD = (Einv * K * D).pruned();
    copyBlock<false>(EinvKD, mat, 0, 0, firstRow<strain>(), firstRow<strain>(),
		     nRows<strain>(), nCols<strain>());

    // We no longer need the D or E matrices, so we can release them
    D = MatrixType();
    E = MatrixType();

    return;
}

void SystemMatrix::addFixedConstraint(int index)
{
    // Add the necessary rows and columns
    addDim<addRow>(mat, subIndex<SubMatrix::fixed>(index) - 1, 1);
    addDim<addCol>(mat, subIndex<SubMatrix::fixed>(index, twistLength) - 1, twistLength);

    return;
}

void SystemMatrix::updateFixedConstraints()
{
    // Start by clearing the relevant submatrices
    clearBlock<SubMatrix::twist, SubMatrix::fixed, true>();
    clearBlock<SubMatrix::strain, SubMatrix::fixed, true>();
    clearBlock<SubMatrix::fixed, SubMatrix::system, true>();

    // Constraint forces are applied to each node with a constraint
    // We add them regardless of whether they are active.
    int index = 0;
    for (auto it = fixedConstraints.cbegin(); it != fixedConstraints.cend(); ++it)
    {
	// Add the force
	int row = subIndex<SubMatrix::fixed>(it->node, twistLength);
	int col = subIndex<SubMatrix::fixed>(index, twistLength);
	for (int i = 0; i < twistLength; i++)
	{
	    mat.coeffRef(row + i, col + i) = 1;
	}
	// If active, add the location
	if (it->active)
	{
	    row = rowIndex<fixedConstraintLocation>(index, twistLength);
	    col = colIndex<fixedConstraintLocation>(it->node, twistLength);
	    for (int i = 0; i < twistLength; i++)
	    {
		mat.coeffRef(row + i, col + i) = 1;
	    }
	}
	// Otherwise, add it to the deactivation matrix
	else
	{
	    row = rowIndex<fixedConstraintActive>(index, twistLength);
	    col = colIndex<fixedConstraintActive>(it->node, twistLength);
	    for (int i = 0; i < twistLength; i++)
	    {
		mat.coeffRef(row + i, col + i) = 1;
	    }
	}
	index++;
    }
    return;
}

template <SystemMatrix::SubMatrix sub, bool prune>
void SystemMatrix::clearBlock()
{
    cossolve::clearBlock<prune>(mat, firstRow<sub>(), firstCol<sub>(), nRows<sub>(), nCols<sub>());
    return;
}

void SystemMatrix::initStiffnessMatrix(Eigen::Ref<const TwistType> diag)
{
    // For now, this just gets put in the EinvKD matrix
    K = MatrixType(nRows<strain>(), nCols<strain>());
    for (int i = firstCol<strain>(); i <= lastCol<strain>(); i++)
    {
	K.coeffRef(i, i) = diag(i % twistLength);
    }
    return;
}

void SystemMatrix::initDerivativeMatrix(const DerivativeMatrixCoefficients& coeffs)
{
    // We make a temporary copy before multiplying it with the stiffness matrix
    D = MatrixType(nRows<strain>(), nCols<strain>());
    int i = firstRow<strain>();
    // First node
    for (; i < rowIndex<strain>(1, twistLength); i++)
    {
	D.coeffRef(i, i) = coeffs.firstNodeCurrentCoeff;
	D.coeffRef(i, i + twistLength) = coeffs.firstNodeNextCoeff;
    }
    // Inner nodes
    for (; i < rowIndex<strain>(-1, twistLength); i++)
    {
	D.coeffRef(i, i - twistLength) = coeffs.innerNodePrevCoeff;
	D.coeffRef(i, i) = coeffs.innerNodeCurrentCoeff;
	D.coeffRef(i, i + twistLength) = coeffs.innerNodeNextCoeff;
    }
    // Last node
    for (; i <= lastRow<strain>(); i++)
    {
	D.coeffRef(i, i - twistLength) = coeffs.lastNodePrevCoeff;
	D.coeffRef(i, i) = coeffs.lastNodeCurrentCoeff;
    }
    return;
}

void SystemMatrix::initIntegralMatrix(const IntegralMatrixCoefficients& coeffs)
{
    // We don't actually add this to the system matrix yet
    E = MatrixType(nRows<strain>(), nCols<strain>());
    int i = firstRow<strain>();
    // First node
    for (; i < rowIndex<strain>(1, twistLength); i++)
    {
	E.coeffRef(i, i) = coeffs.firstNodeCurrentCoeff;
	E.coeffRef(i, i + twistLength) = coeffs.firstNodeNextCoeff;
    }
    // Inner nodes
    for (; i < rowIndex<strain>(-1, twistLength); i++)
    {
	int j = firstCol<strain>() + i % twistLength;
	E.coeffRef(i, j) = coeffs.innerNodeFirstCoeff;
	j += twistLength;
	for (; j < i - twistLength; j += twistLength)
	{
	    E.coeffRef(i, j) = coeffs.innerNodeIntermediateCoeff;
	}
	// Make sure we don't overwrite the first node
	if (i - twistLength >= twistLength)
	{
	    E.coeffRef(i, i - twistLength) = coeffs.innerNodePrevCoeff;
	}
	E.coeffRef(i, i) = coeffs.innerNodeCurrentCoeff;
	E.coeffRef(i, i + twistLength) = coeffs.innerNodeNextCoeff;
    }
    // Last node
    for (; i <= lastRow<strain>(); i++)
    {
	int j = firstCol<strain>() + i % twistLength;
	E.coeffRef(i, j) = coeffs.lastNodeFirstCoeff;
	j += twistLength;
	for (; j < i - twistLength; j += twistLength)
	{
	    E.coeffRef(i, j) = coeffs.lastNodeIntermediateCoeff;
	}
	E.coeffRef(i, i - twistLength) = coeffs.lastNodePrevCoeff;
	E.coeffRef(i, i) = coeffs.lastNodeCurrentCoeff;
    }
    return;
}

} // namespace cossolve
