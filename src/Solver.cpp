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
    logger.log<LogLevel::info>() << "Initializing solver with " << nNodes << " nodes.";
    logger.startTimer("initialization");
    dt = 0.01;
    initSystemMatrices();
    initStateVectors();
    addDistributedForce(0, 1, 9.81*linearDensity*SingleVectorType {0, 0, -1, 0, 0, 0}, false);
    logger.endTimer();

    return;
}

Solver::~Solver()
{
}

void Solver::addForce(ScalarType s, Eigen::Ref<const SingleVectorType> force, bool bodyFrame)
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

void Solver::addPointForce(ScalarType s, Eigen::Ref<const SingleVectorType> force, bool bodyFrame)
{
    // If the force is being placed at the end, we need to double it
    // since the end covers only half a segment.
    int node = std::round(s / ds);
    if (node == 0 || node == (nNodes - 1))
    {
	addForce(s, 2*force*(ds/length));
    }
    else
    {
	addForce(s, force*(ds/length));
    }
    return;
}

void Solver::addDistributedForce(ScalarType s1, ScalarType s2,
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

void Solver::timeStep()
{
    computeInertia();
    // applyContactForces();
    applyForces();
    computeRigidKinematics();
    computeNodalAccelerations();
    applyForces();
    computeStrains();
    computeNodalVelocities();
    solveCoords();

    return;
}

void Solver::initStateVectors()
{
    strains = VectorType::Zero(Sizes::entriesPerVector * nDims);
    freeStrains = VectorType::Zero(Sizes::entriesPerVector * nDims);
    forces = VectorType::Zero(Sizes::entriesPerVector * nDims);
    gBody = std::vector<CoordType>(nNodes);
    gBodyCentroid = std::vector<CoordType>(nNodes);
    gCentroidBody = std::vector<CoordType>(nNodes);
    rigidBodyVelocity = SingleVectorType::Zero();
    rigidBodyAcceleration = SingleVectorType::Zero();
    nodalVelocities = VectorType::Zero(Sizes::entriesPerVector * nDims);
    nodalAccelerations = VectorType::Zero(Sizes::entriesPerVector * nDims);
    
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

    // Compute the initial centroid and moment of inertia
    computeInertia();

    return;
}

void Solver::initSystemMatrices()
{
    initStiffnessMatrix();
    initInertiaMatrix();
    initDerivativeMatrix();
    initIntegralMatrix();
    initStrainEqnMatrix();
    
    // Create adjoint matrices (we don't need to fill it since it needs to be regenerated at each
    // solution step anyways)
    Af = MatrixType(Sizes::entriesPerVector * nDims, Sizes::entriesPerVector * nDims);
    Av = MatrixType(Af.rows(), Af.cols());
    AfK = MatrixType(Af.rows(), Af.cols());
    AvM = MatrixType(Af.rows(), Af.cols());

    return;
}

void Solver::initStiffnessMatrix()
{
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
    for (i = 0; i < nodeIndex(-1); i++)
    {
	tripletList.emplace_back(i, i, Kdiag(i % Sizes::entriesPerVector));
    }
    // Bottom right block is identity (for BCs)
    for (; i < endIndex; i++)
    {
	tripletList.emplace_back(i, i, 1);
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
    return;
}

void Solver::initInertiaMatrix()
{
    // Assuming a rectangular cross section, we calculate the moments of inertia
    ScalarType height = sqrt(12*momentY/area);
    ScalarType width = area / height;
    SingleVectorType Mdiag = SingleVectorType::Zero();
    Mdiag(0) = linearDensity;
    Mdiag(1) = Mdiag(0);
    Mdiag(2) = Mdiag(1);
    Mdiag(3) = linearDensity/area * momentX;
    Mdiag(4) = linearDensity/12 * (height*height + length*length / (nSegments*nSegments));
    Mdiag(5) = linearDensity/12 * (width*width + length*length / (nSegments*nSegments));

    M = MatrixType(Sizes::entriesPerVector * nDims, Sizes::entriesPerVector * nDims);
    std::vector<Eigen::Triplet<ScalarType>> tripletList;
    tripletList.reserve(Sizes::entriesPerVector * nDims);
    int i;
    // The first and last block get half values
    for (i = 0; i < nodeIndex(1); i++)
    {
	tripletList.emplace_back(i, i, Mdiag(i) / 2);
    }
    for (; i < nodeIndex(-2); i++)
    {
	tripletList.emplace_back(i, i, Mdiag(i % Sizes::entriesPerVector));
    }
    for (; i < nodeIndex(-1); i++)
    {
	tripletList.emplace_back(i, i, Mdiag(i % Sizes::entriesPerVector) / 2);
    }
    M.setFromTriplets(tripletList.begin(), tripletList.end());

    // Also initialize the rigid body inertia tensor
    Mrigid.setZero();
    return;
}

void Solver::initDerivativeMatrix()
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

void Solver::initIntegralMatrix()
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

void Solver::initStrainEqnMatrix()
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

void Solver::generateAdjointMatrix()
{
    SingleMatrixType adf;
    SingleMatrixType adv;
    for (int i = 0; i < nNodes; i++)
    {
	adf.setZero();
	adv.setZero();
	adjoint(nodeVectorConst(strains, i), adf);
	adjoint(nodeVectorConst(nodalVelocities, i), adv);
	adf.transposeInPlace();
	adv.transposeInPlace();
	for (int j = 0; j < Sizes::entriesPerVector; j++)
	{
	    for (int k = 0; k < Sizes::entriesPerVector; k++)
	    {
		Af.coeffRef(nodeIndex(i)+j, nodeIndex(i)+k) = adf(j, k);
		Av.coeffRef(nodeIndex(i)+j, nodeIndex(i)+k) = adv(j, k);
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
	    g.block<3, 3>(0, 0) = gBody[node].block<3, 3>(0, 0);
	    Adg.setZero();
	    Adjoint(g, Adg);
	    nodeVector(forces, node) = Adg.transpose()*force;
	}
    }

    // The last vector in our forces vector holds the boundary condition (strain at the imaginary (n+1)th node)
    SingleMatrixType invK = SingleMatrixType(K.topLeftCorner(Sizes::entriesPerVector,
							     Sizes::entriesPerVector)).inverse();
    SingleMatrixType Mblock = M.block(0, 0, Sizes::entriesPerVector, Sizes::entriesPerVector);
    nodeVector(forces, -1) = -(nodeVector(strains, -2) - ds/2*invK*((nodeVector(forces, -2)) - Mblock*nodeVector(nodalAccelerations, -2)));

    return; 
}

void Solver::applyContactForces()
{
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

void Solver::computeStrains()
{
    logger.startTimer("Solve Strain");
    
    // Reconstruct the adjoint matrices
    logger.startTimer("Adjoint Matrix");
    generateAdjointMatrix();
    AfK = (Af * K).pruned();
    AvM = (Av * M).pruned();
    logger.endTimer();
    
    // Compute
    logger.startTimer("compute");
    sparseLuSolver.compute((EinvKD - AfK).pruned());
    logger.endTimer();
    logger.startTimer("compute RHS");
    VectorType rhs = -(AfK * freeStrains) - forces + (M * nodalAccelerations) - (AvM * nodalVelocities);
    logger.endTimer();
    strains = sparseLuSolver.solve(rhs);
    logger.endTimer();
    
    return;
}

void Solver::solveCoords(int baseNode)
{
    logger.startTimer("coordinate calculation");
    // First, update baseNode position due to rigid body velocity
    CoordType vHat = CoordType::Zero();
    se3hat(nodeVectorConst(nodalVelocities, baseNode), vHat);
    gBody[baseNode] *= (vHat * dt).exp();

    // Now apply deformation due to strains
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
    logger.endTimer();
    return;
}

void Solver::computeInertia()
{
    // Assume the mass is located at the midpoint of each segments
    gCentroid.setIdentity();
    CoordType g;
    SingleMatrixType Adg = SingleMatrixType::Zero();
    for (int seg = 0; seg < nSegments; seg++)
    {
	int nearestNode = intermediateTransform((seg + 0.5)*ds, g);
	g = gBody[nearestNode]*g;
	origin(gCentroid) += origin(g);
    }
    origin(gCentroid) /= nSegments;

    // Recompute the transformations between the centroid and body frames
    CoordType gCentroidInv = gCentroid.inverse();
    for (int node = 0; node < nNodes; node++)
    {
	gBodyCentroid[node] = gBody[node].inverse() * gCentroid;
	gCentroidBody[node] = gCentroidInv * gBody[node];
    }
    
    // Moment of inertia is computed by transforming the inertia tensor
    // of each node to the centroid frame and adding them up
    Mrigid.setZero();
    SingleMatrixType Mblock = M.block(nodeIndex(1), nodeIndex(1), Mrigid.rows(), Mrigid.cols());
    auto calcNode = [&](int node)
	{
	    Adg.setZero();
	    Adjoint(gBodyCentroid[node], Adg);
	    return Adg.transpose() * Mblock * Adg;
	};
    Mrigid += calcNode(0) / 2;
    for (int i = 1; i < (nNodes - 1); i++)
    {
	Mrigid += calcNode(i);
    }
    Mrigid += calcNode(nNodes - 1) / 2;
    // The body inertias are per unit length. We need to convert them back to
    // discrete values.
    Mrigid *= length / nSegments;

    return;
}

void Solver::computeRigidKinematics()
{
    // First, calculate the net force by transforming all of the nodal
    // forces to the centroid
    SingleVectorType netForce = SingleVectorType::Zero();
    SingleMatrixType Adg = SingleMatrixType::Zero();
    for (int i = 0; i < nNodes; i++)
    {
	Adjoint(gBodyCentroid[i], Adg);
	netForce += Adg.transpose() * nodeVector(forces, i);
    }
    // The body forces are per unit length. We need to convert back to
    // discrete values.
    netForce *= length / nSegments;

    // Now compute acceleration
    Adg.setZero();
    adjoint(rigidBodyVelocity, Adg);
    rigidBodyAcceleration = Mrigid.inverse() * (netForce + Adg.transpose() * Mrigid * rigidBodyVelocity);

    // Update velocity
    rigidBodyVelocity += rigidBodyAcceleration * dt;
    
    return;
}

void Solver::computeNodalVelocities()
{
    SingleMatrixType Adg = SingleMatrixType::Zero();
    for (int node = 0; node < nNodes; node++)
    {
	Adjoint(gBodyCentroid[node], Adg);
	nodeVector(nodalVelocities, node) = Adg * rigidBodyVelocity;
    }
    return;
}

void Solver::computeNodalAccelerations()
{
    SingleMatrixType Adg = SingleMatrixType::Zero();
    for (int node = 0; node < nNodes; node++)
    {
	Adjoint(gBodyCentroid[node], Adg);
	nodeVector(nodalAccelerations, node) = Adg * rigidBodyAcceleration;
    }
    return;
}

} // namespace cossolve
