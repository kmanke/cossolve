/* Implementation of Solver class.
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

#include "Solver.h"
#include "se3.h"

#include <chrono>
#include <iostream>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

namespace cossolve {

Solver::Solver(int _nNodes, ScalarType _tensileModulus, ScalarType _poissonRatio, ScalarType _shearModulus,
	       ScalarType _momentX, ScalarType _momentY, ScalarType _momentZ, ScalarType _area,
	       ScalarType _length, ScalarType _linearDensity)
    : tensileModulus(_tensileModulus), poissonRatio(_poissonRatio),
      shearModulus(_shearModulus == 0 ? tensileModulus / (2*(1 + poissonRatio)) : _shearModulus),
      momentX(_momentX), momentY(_momentY), momentZ(_momentZ),
      area(_area), length(_length), linearDensity(_linearDensity),
      nNodes(_nNodes), nSegments(nNodes - 1), nDims(nNodes + 1), ds(1.0 / nSegments),
      endIndex(nDims * Sizes::entriesPerVector), logger(std::cout)
{
    logger(Logger::info) << "Initializing solver with " << nNodes << " nodes.";
    logger.startTimer("initialization");
    initSystemMatrices();
    initStateVectors();
    logger.endTimer();
}

Solver::~Solver()
{
}

void Solver::addForce(ScalarType s, Eigen::Ref<const SingleVectorType> force, bool bodyFrame)
{
    // First, transform the force to its nearest node.
    CoordType g;
    ForceEntry newForce;
    g.setZero();
    int node = intermediateTransform(s, g);

    // Calculate the adjoint of the transformation
    SingleMatrixType Adg;
    Adg.setZero();
    Adjoint(g.inverse(), Adg);

    // Fill in the tuple
    std::get<0>(newForce) = Adg.transpose() * force;
    std::get<1>(newForce) = node;
    std::get<2>(newForce) = bodyFrame;

    // Adjust the force from discrete to distributed
    std::get<0>(newForce) /= (ds * length);

    // If the force is at the tip or the base, we need to ScalarType it to account for the
    // piecewise linear derivative behaviour near the ends
    if (node == 0 || node == (nNodes - 1))
    {
	std::get<0>(newForce) *= 2;
    }
    externalForces.emplace_back(std::move(newForce));
    
    return;
}

void Solver::initStateVectors()
{
    strains = VectorType::Zero(Sizes::entriesPerVector * nDims);
    freeStrains = VectorType::Zero(Sizes::entriesPerVector * nDims);
    forces = VectorType::Zero(Sizes::entriesPerVector * nDims);
    gs = std::vector<CoordType>(nNodes);

    Eigen::Vector<ScalarType, 4> p {0, 0, 0, 1};
    for (auto gsIt = gs.begin(); gsIt != gs.end(); ++gsIt)
    {
	*gsIt = CoordType::Zero();
	(*gsIt).block<3, 3>(0, 0).setIdentity();
	(*gsIt).block<4, 1>(0, 3) = p;

	p(0) += length / nSegments;
    }

    // Compute free strains based on initial geometry
    for (int i = 0; i < (nNodes - 1); i++)
    {
	CoordType dg = (gs[i+1] - gs[i]) / ds;
	se3unhat(gs[i].inverse()*dg, nodeVector(freeStrains, i));
    }
    // Set the free strain of the end to be the same as the previous node
    nodeVector(freeStrains, -2) = nodeVector(freeStrains, -3);
    nodeVector(freeStrains, -1) = nodeVector(freeStrains, -3);
    
    // Initial strains are just set to free strains
    strains = freeStrains;

    // Compute the initial centroid and moment of inertia
    computeInertia();

    return;
}

void Solver::initSystemMatrices()
{
    // Prepare the stiffness matrix
    SingleVectorType Kdiag = SingleVectorType::Zero();
    Kdiag(0) = tensileModulus * area;
    Kdiag(1) = shearModulus * area;
    Kdiag(2) = Kdiag(1);
    Kdiag(3) = shearModulus * momentX;
    Kdiag(4) = tensileModulus * momentY;
    Kdiag(5) = tensileModulus * momentZ;
    
    K = MatrixType(Sizes::entriesPerVector*nDims, Sizes::entriesPerVector*nDims);
    std::vector<Eigen::Triplet<ScalarType>> tripletList;
    tripletList.reserve(Sizes::entriesPerVector * nDims);
    int i;
    for (i = 0; i < nNodes; i++)
    {
	for (int j = 0; j < Kdiag.size(); j++)
	{
	    tripletList.emplace_back(nodeIndex(i) + j, nodeIndex(i) + j,
				     Kdiag(j));
	}
    }
    // Bottom right block is identity (for BCs)
    for (i = nodeIndex(-1); i < endIndex; i++)
    {
	tripletList.emplace_back(i, i, 1);
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());

    // Prepare the derivative matrix
    D = MatrixType(Sizes::entriesPerVector*nDims, Sizes::entriesPerVector*nDims);
    tripletList.clear();
    tripletList.reserve(Sizes::entriesPerVector*2*nDims);
    // The first block is special
    i = 0;
    for (; i < nodeIndex(1); i++)
    {
	tripletList.emplace_back(i, i, -1/ds);
	tripletList.emplace_back(i, i+Sizes::entriesPerVector, 1/ds);
    }
    // Interior blocks
    for (; i < nodeIndex(-1); i++)
    {
	tripletList.emplace_back(i, i-Sizes::entriesPerVector, -0.5/ds);
	tripletList.emplace_back(i, i+Sizes::entriesPerVector, 0.5/ds);
    }
    // Bottom 6 rows are special for BCs
    for (; i < endIndex; i++)
    {
	tripletList.emplace_back(i, i, 1);
    }
    D.setFromTriplets(tripletList.begin(), tripletList.end());

    // Prepare the integral matrix
    E = MatrixType(Sizes::entriesPerVector*nDims, Sizes::entriesPerVector*nDims);
    tripletList.clear();
    tripletList.reserve(Sizes::entriesPerVector*3*nDims);
    // The first 6 rows are special
    i = 0;
    for (; i < nodeIndex(1); i++)
    {
	tripletList.emplace_back(i, i, 0.5);
	tripletList.emplace_back(i, i+Sizes::entriesPerVector, 0.5);
    }
    // Interior blocks
    for (; i < nodeIndex(-2); i++)
    {
	tripletList.emplace_back(i, i-Sizes::entriesPerVector, 0.25);
	tripletList.emplace_back(i, i, 0.5);
	tripletList.emplace_back(i, i+Sizes::entriesPerVector, 0.25);
    }
    // Bottom 12 rows are special
    for (; i < nodeIndex(-1); i++)
    {
	tripletList.emplace_back(i, i-Sizes::entriesPerVector, 0.25);
	tripletList.emplace_back(i, i, 0.5);
    }
    for (; i < endIndex; i++)
    {
	tripletList.emplace_back(i, i, 1);
    }
    E.setFromTriplets(tripletList.begin(), tripletList.end());

    // Generate KD (used to generate EinvKD)
    MatrixType KD = (K*D).pruned();
    /*for (i = nodeIndex(-2); i < nodeIndex(-1); i++)
    {
	KD.coeffRef(i, i) = -1;
	KD.coeffRef(i, i + Sizes::entriesPerVector) = 1;
    }
    for (; i < endIndex; i++)
    {
	KD.coeffRef(i, i) = 1;
	}*/
    // Now generate EinvKD
    EinvKD = MatrixType(Sizes::entriesPerVector*nDims, Sizes::entriesPerVector*nDims);
    MatrixType I(EinvKD.rows(), EinvKD.cols());
    I.setIdentity();
    matrixSolver.compute(E);
    EinvKD = (matrixSolver.solve(I) * KD).pruned(1, 1e-2);

    // Create A and AK (we don't need to fill it since it needs to be regenerated at each
    // solution step anyways)
    A = MatrixType(Sizes::entriesPerVector * nDims, Sizes::entriesPerVector * nDims);
    AK = MatrixType(A.rows(), A.cols());
    
    return;
}

void Solver::generateAdjointMatrix()
{
    SingleMatrixType adf = SingleMatrixType::Zero();
    for (int i = 0; i < nNodes; i++)
    {
	adf.setZero();
	adjoint(nodeVectorConst(strains, i), adf);
	adf.transposeInPlace();
	for (int j = 0; j < Sizes::entriesPerVector; j++)
	{
	    for (int k = 0; k < Sizes::entriesPerVector; k++)
	    {
		A.coeffRef(nodeIndex(i)+j, nodeIndex(i)+k) = adf(j, k);
	    }
	}
    }
    return;
}

/*void Solver::setupBoundaryCondition()
{
    // The boundary condition is represented by the last two equations in the matrix
    // which basically state:
    // (1) (fn+1 - fn-1)/ds = Kf'n (calculated)
    // (2) fn+1 = fn + ds*Kf'n (constant)
    // This essentially just creates an imaginary point just past the end of the
    // matrix 
    }*/

void Solver::applyForces()
{
    CoordType g;
    SingleMatrixType Adg;

    // First apply external forces
    for (auto it = externalForces.cbegin(); it != externalForces.cend(); ++it)
    {
	Eigen::Ref<const SingleVectorType> force = std::get<0>(*it);
	const int& node = std::get<1>(*it);
	const bool& bodyFrame = std::get<2>(*it);
	
	if (bodyFrame) // Force can be applied directly with no transformation
	{
	    nodeVector(forces, node) = force;
	}
	else // Need to rotate from the spatial frame into the body frame
	{
	    // Construct a coordinate transformation which only rotates
	    g.setIdentity();
	    g.block<3, 3>(0, 0) = gs[node].block<3, 3>(0, 0);
	    Adg.setZero();
	    Adjoint(g, Adg);
	    nodeVector(forces, node) = Adg.transpose()*force;
	}
    }

    // The last vector in our forces vector holds the boundary condition (strain at the imaginary (n+1)th node)
    SingleMatrixType invK = SingleMatrixType(K.topLeftCorner(Sizes::entriesPerVector,
							     Sizes::entriesPerVector)).inverse();
    nodeVector(forces, -1) = -(nodeVector(strains, -2) - ds/2*invK*nodeVector(forces, -2));

    return; 
}

int Solver::intermediateTransform(ScalarType s, Eigen::Ref<CoordType> g) const
{
    // First, find the nearest node to s
    int nearestNode = std::round(s / ds);
    ScalarType sk = nearestNode * ds;
    
    // Compute the strain integral
    SingleVectorType integratedStrain = SingleVectorType::Zero();
    if (s > sk)
    {
	intermediateStrainIntegral(sk, s, integratedStrain);
    }
    else
    {
	intermediateStrainIntegral(s, sk, integratedStrain);
	integratedStrain = -integratedStrain;
    }
    CoordType integratedStrainHat;
    integratedStrainHat.setZero();
    se3hat(integratedStrain, integratedStrainHat);

    // Compute the result
    g = integratedStrainHat.exp();
    return nearestNode;
}

void Solver::intermediateStrainIntegral(ScalarType s1, ScalarType s2, Eigen::Ref<SingleVectorType> result) const
{
    // Calculate the first node before s1 and the first node before s2
    int node1 = std::floor(s1 / ds);
    int node2 = std::floor(s2 / ds);

    // First, integrate from node1 to node2
    nodalStrainIntegral(node1, node2, result);

    // Now we need to fix up the start and the end
    relativeStrainIntegral(node1, s1 - node1*ds, result, true); // Subtract the first part
    relativeStrainIntegral(node2, s2 - node2*ds, result, false); // Add the end part

    return;
}

void Solver::nodalStrainIntegral(int node1, int node2, Eigen::Ref<SingleVectorType> result) const
{
    while (node1 < node2)
    {
	if (node1 == 0) // First node
	{
	    Eigen::Ref<const SingleVectorType> f0 = nodeVectorConst(strains, 0);//strains.block<6, 1>(0, 0);
	    Eigen::Ref<const SingleVectorType> f1 = nodeVectorConst(strains, 1);
	    Eigen::Ref<const SingleVectorType> f2 = nodeVectorConst(strains, 2);
	    result += ds*(5.0/12.0*f0 + 8.0/12.0*f1 - 1.0/12.0*f2);
	}
	else // General case
	{
	    Eigen::Ref<const SingleVectorType> fprev = nodeVectorConst(strains, node1 - 1);
	    Eigen::Ref<const SingleVectorType> fk = nodeVectorConst(strains, node1);
	    Eigen::Ref<const SingleVectorType> fnext = nodeVectorConst(strains, node1 + 1);
	    Eigen::Ref<const SingleVectorType> fnext2 = nodeVectorConst(strains, node1 + 2);
	    result += ds*(1.0/6.0*(fnext - fprev) + 11.0/12.0*fk + 1.0/12.0*fnext2);
	}
	node1++;
    }
    return;
}

void Solver::relativeStrainIntegral(int node, ScalarType x, Eigen::Ref<SingleVectorType> result, bool subtract) const
{
    if (node == 0) // Special case for the first segment
    {
	Eigen::Ref<const SingleVectorType> f0 = nodeVectorConst(strains, 0);
	Eigen::Ref<const SingleVectorType> f1 = nodeVectorConst(strains, 1);
	Eigen::Ref<const SingleVectorType> f2 = nodeVectorConst(strains, 2);

	result += (subtract ? -1 : 1) * (x*f0 + x*x/(2*ds)*(f1 - f0) + x*x*x/(12*ds*ds)*(f2 + 2*f1 - 3*f0));
    }
    else if (node < (nNodes - 1)) // General case for middle segments
    {
	Eigen::Ref<const SingleVectorType> fprev = nodeVectorConst(strains, node - 1);
	Eigen::Ref<const SingleVectorType> fk = nodeVectorConst(strains, node);
	Eigen::Ref<const SingleVectorType> fnext = nodeVectorConst(strains, node + 1);
	Eigen::Ref<const SingleVectorType> fnext2 = nodeVectorConst(strains, node + 2);
	result += (subtract ? -1 : 1) * (x*fk + x*x/(4*ds)*(fnext - fprev) + x*x*x/(12*ds*ds)*(fnext2 - fk + fnext - fprev));
    }
    // If node >= nNodes - 1, then we're at the tip, so do nothing at all
    return;
}

void Solver::solveStrains()
{
    logger.startTimer("Solve Strain");
    // First, update the applied forces
    applyForces();

    // Reconstruct the adjoint matrix and recalculate AK
    generateAdjointMatrix();
    AK = (A * K).pruned();
    
    // Compute
    strainSolver.compute((EinvKD - AK).pruned());
    VectorType rhs = -(AK * freeStrains) - forces;
    strains = strainSolver.solve(rhs);
    logger.endTimer();
    
    return;
}

void Solver::solveCoords()
{
    logger.startTimer("coordinate calculation");
    SingleVectorType integratedStrain;
    CoordType integratedStrainHat = CoordType::Zero();
    for (int i = 1; i < nNodes; i++)
    {
	integratedStrain.setZero();
	nodalStrainIntegral(i - 1, i, integratedStrain);
	se3hat(integratedStrain, integratedStrainHat);
	gs[i] = gs[i-1]*integratedStrainHat.exp();
    }
    computeInertia();
    logger.endTimer();
    return;
}

void Solver::computeInertia()
{
    // Option 1: Mass located at midpoint of each segment
    centroid.setZero();
    CoordType segCentre;
    for (int seg = 0; seg < nSegments; seg++)
    {
	int nearestNode = intermediateTransform((seg + 0.5)*ds, segCentre);
	segCentre = gs[nearestNode]*segCentre;
	origin(centroid) += origin(segCentre);
    }
    origin(centroid) /= nSegments;

    return;
}

} // namespace cossolve
