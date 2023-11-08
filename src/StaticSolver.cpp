/* Implementation of StaticSolver class.
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

#include "StaticSolver.h"
#include "se3.h"

#include <chrono>
#include <iostream>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

namespace cossolve {

StaticSolver::StaticSolver(int _nNodes, ScalarType _tensileModulus, ScalarType _poissonRatio, ScalarType _shearModulus,
			   ScalarType _momentX, ScalarType _momentY, ScalarType _momentZ, ScalarType _area,
			   ScalarType _length, ScalarType _linearDensity)
    : tensileModulus(_tensileModulus), poissonRatio(_poissonRatio),
      shearModulus(_shearModulus == 0 ? tensileModulus / (2*(1 + poissonRatio)) : _shearModulus),
      momentX(_momentX), momentY(_momentY), momentZ(_momentZ),
      area(_area), length(_length), linearDensity(_linearDensity),
      nNodes(_nNodes), nSegments(nNodes - 1), nDims(nNodes + 1), ds(1.0 / nSegments),
      endIndex(nDims * Sizes::entriesPerVector), logger(std::cout)
{
    logger.log<LogLevel::info>() << "Initializing solver with " << nNodes << " nodes.";
    logger.startTimer("initialization");
    dt = 0.00001;
    t = 0;
    initSystemMatrices();
    initStateVectors();
//    addDistributedForce(0, 1, 9.81*linearDensity*SingleVectorType {0, 0, -1, 0, 0, 0}, false);

    logger.endTimer();

    return;
}

StaticSolver::~StaticSolver()
{
}

void StaticSolver::addForce(ScalarType s, Eigen::Ref<const SingleVectorType> force, bool bodyFrame)
{
    // First, transform the force to its nearest node
    CoordType g;
    ForceEntry newForce;
    g.setZero();
    int node = intermediateTransform(s, g);

    // Calculate the adjoint of the transformation
    SingleMatrixType Adg;
    Adg.setZero();
    Adjoint(g.inverse(), Adg);

    // Fill in the tuple and move it over
    std::get<0>(newForce) = Adg.transpose() * force;
    std::get<1>(newForce) = node;
    std::get<2>(newForce) = bodyFrame;
    logger.log() << std::get<0>(newForce);
    externalForces.emplace_back(std::move(newForce));
    
    return;
}

void StaticSolver::addPointForce(ScalarType s, Eigen::Ref<const SingleVectorType> force, bool bodyFrame)
{
    // If the force is being placed at the end, we need to double it
    // since the end covers only half a segment.
    int node = std::round(s / ds);
    if (node == 0 || node == (nNodes - 1))
    {
	addForce(s, 2*force/(ds*length));
    }
    else
    {
	addForce(s, force/(ds*length));
    }
    return;
}

void StaticSolver::addDistributedForce(ScalarType s1, ScalarType s2,
				 Eigen::Ref<const SingleVectorType> force, bool bodyFrame)
{
    // First deal with the leftmost node
    int leftNode = std::round(s1 / ds);
    ScalarType sNode = leftNode * ds;
    ScalarType width = sNode + ds/2 - s1;
    ScalarType mid = s1 + width / 2;
    addForce(mid, force * width / ds, bodyFrame);
    //addForce(sNode, force * width / ds, bodyFrame);
    
    // Now the rightmost node
    int rightNode = std::round(s2 / ds);
    sNode = rightNode * ds;
    width = sNode + ds/2 - s2;
    mid = s2 - width / 2;
    addForce(mid, force * width / ds, bodyFrame);
    //addForce(sNode, force * width / ds, bodyFrame);
    
    // Interior nodes
    for (int node = leftNode + 1; node < rightNode; node++)
    {
	addForce(ds * node, force, bodyFrame);
    }

    return;
}

void StaticSolver::timeStep()
{
    applyForces(1);
    applyContactForces();
    computeStrains();
    computeDeformation();

    return;
}

void StaticSolver::initStateVectors()
{
    strains = VectorType::Zero(Sizes::entriesPerVector * nDims);
    freeStrains = VectorType::Zero(Sizes::entriesPerVector * nDims);
    forces = VectorType::Zero(Sizes::entriesPerVector * nDims);
    gBody = std::vector<CoordType>(nNodes);
    contactNodes.reserve(nNodes);
    
    Eigen::Vector<ScalarType, 4> p {0, 0, 0, 1};
    for (auto gBodyIt = gBody.begin(); gBodyIt != gBody.end(); ++gBodyIt)
    {
	*gBodyIt = CoordType::Zero();
	(*gBodyIt).block<3, 3>(0, 0).setIdentity();
	(*gBodyIt).block<4, 1>(0, 3) = p;

	p(0) += length / nSegments;
    }

    // Compute free strains based on initial geometry
    for (int i = 0; i < (nNodes - 1); i++)
    {
	CoordType dg = (gBody[i+1] - gBody[i]) / ds;
	se3unhat(gBody[i].inverse()*dg, nodeVector(freeStrains, i));
    }
    // Set the free strain of the end to be the same as the previous node
    nodeVector(freeStrains, -2) = nodeVector(freeStrains, -3);
    nodeVector(freeStrains, -1) = nodeVector(freeStrains, -3);
    
    // Initial strains are just set to free strains
    strains = freeStrains;

    return;
}

void StaticSolver::initSystemMatrices()
{
    initStiffnessMatrix();
    initDerivativeMatrix();
    initIntegralMatrix();
    initStrainEqnMatrix();
    
    // Create adjoint matrices (we don't need to fill it since it needs to be regenerated at each
    // solution step anyways)
    Af = MatrixType(Sizes::entriesPerVector * nDims, Sizes::entriesPerVector * nDims);
    AfK = MatrixType(Af.rows(), Af.cols());

    return;
}

void StaticSolver::initStiffnessMatrix()
{
    SingleVectorType Kdiag = SingleVectorType::Zero();
    Kdiag(0) = tensileModulus * area;
    Kdiag(1) = shearModulus * area;
    Kdiag(2) = Kdiag(1);
    Kdiag(3) = shearModulus * momentX;
    Kdiag(4) = tensileModulus * momentY;
    Kdiag(5) = tensileModulus * momentZ;

    K = MatrixType(Sizes::entriesPerVector*nNodes, Sizes::entriesPerVector*nNodes);
    std::vector<Eigen::Triplet<ScalarType>> tripletList;
    tripletList.reserve(Sizes::entriesPerVector * nNodes);
    int i;
    for (i = 0; i < endIndex; i++)
    {
	tripletList.emplace_back(i, i, Kdiag(i % Sizes::entriesPerVector));
    }
    return;
}

void StaticSolver::initDerivativeMatrix()
{
    D = MatrixType(Sizes::entriesPerVector*nDims, Sizes::entriesPerVector*nDims);
    std::vector<Eigen::Triplet<ScalarType>> tripletList;
    tripletList.reserve(Sizes::entriesPerVector*nDims*2);
    // The first block is special
    int i;
    for (i = 0; i < nodeIndex(1); i++)
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
    return;
}

void StaticSolver::initIntegralMatrix()
{
    E = MatrixType(Sizes::entriesPerVector*nDims, Sizes::entriesPerVector*nDims);
    std::vector<Eigen::Triplet<ScalarType>> tripletList;
    tripletList.reserve(Sizes::entriesPerVector*3*nDims);
    // The first 6 rows are special
    int i;
    for (i = 0; i < nodeIndex(1); i++)
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
    return;
}

void StaticSolver::initStrainEqnMatrix()
{
    // Generate KD (used to generate EinvKD)
    MatrixType KD = (K*D).pruned();
    EinvKD = MatrixType(Sizes::entriesPerVector*nDims, Sizes::entriesPerVector*nDims);
    MatrixType I(EinvKD.rows(), EinvKD.cols());
    I.setIdentity();
    sparseLuSolver.compute(E);
    EinvKD = (sparseLuSolver.solve(I) * KD).pruned(1, 1e-2);

    return;
}

void StaticSolver::generateAdjointMatrix()
{
    SingleMatrixType adf;
    for (int i = 0; i < nNodes; i++)
    {
	adf.setZero();
	adjoint(nodeVectorConst(strains, i), adf);
	adf.transposeInPlace();
	for (int j = 0; j < Sizes::entriesPerVector; j++)
	{
	    for (int k = 0; k < Sizes::entriesPerVector; k++)
	    {
		Af.coeffRef(nodeIndex(i)+j, nodeIndex(i)+k) = adf(j, k);
	    }
	}
    }
    return;
}

void StaticSolver::applyForces(ScalarType scaleFactor)
{
    CoordType g;
    SingleMatrixType Adg;
    forces.setZero();
    
    // Apply all external forces
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
	    g.block<3, 3>(0, 0) = gBody[node].block<3, 3>(0, 0);
	    Adg.setZero();
	    Adjoint(g, Adg);
	    nodeVector(forces, node) = Adg.transpose() * force * scaleFactor;
	}
    }

    return; 
}

void StaticSolver::applyContactForces()
{
    static constexpr ScalarType floorHeight = -100;
    // Find areas which are below the ground
    contactNodes.clear();
    for (int node = 0; node < nNodes; node++)
    {
	ScalarType clearance = origin(gBody[node])(2) - floorHeight;
	if (clearance < 0)
	{
	    logger.log() << "Contact detected at node " << node;
	    contactNodes.emplace_back(node);
	    // Shift the point until it is coincident with the floor
	    origin(gBody[node])(2) = floorHeight;
	}
    }
    VectorType impliedStrains = VectorType::Zero(Sizes::entriesPerVector * nDims);
    computeStrainFromCoords(impliedStrains);
    nodeVector(impliedStrains, -1) = nodeVectorConst(strains, -1);
    nodeVector(impliedStrains, -2) = nodeVectorConst(strains, -2);
    // Recalculate forces implied by these strains
    //VectorType impliedForces = -((EinvKD - AfK)*impliedStrains + AfK*freeStrains);
    VectorType impliedForces = VectorType::Zero(Sizes::entriesPerVector * nDims);
    for (int node = 0; node < nNodes; node++)
    {
	nodeVector(impliedForces, node) = AfK.block(nodeIndex(node), nodeIndex(node), 6, 6)*(nodeVector(impliedStrains, node) - nodeVector(freeStrains, node))
	    - K.block(nodeIndex(node), nodeIndex(node), 6, 6)*(nodeVector(impliedStrains, node+1) - nodeVector(impliedStrains, node))/ds;
    }
    // For each node in contact, replace the total force with the implied force
    for (auto it = contactNodes.cbegin(); it != contactNodes.cend(); ++it)
    {
	nodeVector(forces, *it) = nodeVectorConst(impliedForces, *it);
    }
    Eigen::MatrixX<ScalarType> comparison = Eigen::MatrixX<ScalarType>(6 * nDims, 4);
    comparison.col(0) = strains;
    comparison.col(1) = impliedStrains;
    comparison.col(2) = strains - impliedStrains;
    comparison.col(3) = impliedForces;
    logger.log() << "Implied strains:\n" << comparison;

    return;
}

int StaticSolver::intermediateTransform(ScalarType s, Eigen::Ref<CoordType> g) const
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

void StaticSolver::intermediateStrainIntegral(ScalarType s1, ScalarType s2, Eigen::Ref<SingleVectorType> result) const
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

void StaticSolver::nodalStrainIntegral(int node1, int node2, Eigen::Ref<SingleVectorType> result) const
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

void StaticSolver::relativeStrainIntegral(int node, ScalarType x, Eigen::Ref<SingleVectorType> result, bool subtract) const
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

void StaticSolver::computeStrains()
{
    logger.startTimer("Solve Strain");
    // Reconstruct the adjoint matrices
    logger.startTimer("Adjoint Matrix");
    generateAdjointMatrix();
    AfK = (Af * K).pruned();
    logger.endTimer();

    // The last vector in our forces vector holds the boundary condition (strain at the imaginary (n+1)th node)
    SingleMatrixType invK = SingleMatrixType(K.topLeftCorner(Sizes::entriesPerVector,
							     Sizes::entriesPerVector)).inverse();
    nodeVector(forces, -1) = -(nodeVector(strains, -2) - ds/2*invK*(nodeVector(forces, -2)));
    
    // Compute
    logger.startTimer("compute");
    sparseLuSolver.compute((EinvKD - AfK).pruned());
    logger.endTimer();
    logger.startTimer("compute RHS");
    VectorType rhs = -(AfK * freeStrains) - forces;
    logger.endTimer();
    strains = sparseLuSolver.solve(rhs);
    logger.endTimer();
    
    return;
}

void StaticSolver::computeDeformation(int baseNode)
{
    SingleVectorType integratedStrain;
    CoordType integratedStrainHat = CoordType::Zero();
    for (int node = baseNode - 1; node >= 0; node--)
    {
	integratedStrain.setZero();
	nodalStrainIntegral(node, node + 1, integratedStrain);
	se3hat(integratedStrain, integratedStrainHat);
	gBody[node] = gBody[node+1]*(-integratedStrainHat.exp());
    }
    for (int node = baseNode + 1; node < nNodes; node++)
    {
	integratedStrain.setZero();
	nodalStrainIntegral(node - 1, node, integratedStrain);
	se3hat(integratedStrain, integratedStrainHat);
	gBody[node] = gBody[node-1]*integratedStrainHat.exp();
    }
    return;
}

void StaticSolver::computeStrainFromCoords(Eigen::Ref<VectorType> dest)
{
    // Do a forward difference for the first
    CoordType strainHat = CoordType::Zero();
    strainHat = gBody[0].inverse() * (gBody[1] - gBody[0]) / ds;
    se3unhat(strainHat, nodeVector(dest, 0));
    for (int node = 1; node < (nNodes - 1); node++)
    {
	strainHat = gBody[node].inverse() * (gBody[node+1] - gBody[node-1]) / (2 * ds);
	se3unhat(strainHat, nodeVector(dest, node));
    }
    strainHat = gBody[nNodes-1].inverse() * (gBody[nNodes-1] - gBody[nNodes-2]) / ds;
    se3unhat(strainHat, nodeVector(dest, nNodes-1));

    return;
}

} // namespace cossolve
