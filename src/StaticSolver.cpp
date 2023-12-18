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
#include "LieGroupOps.h"
#include "MatrixOps.h"
#include "Indexer.h"
#include "Precomputed.h"

#include <eigen3/unsupported/Eigen/MatrixFunctions>

#include <iostream>

namespace cossolve {

StaticSolver::StaticSolver(SolverParameters&& params)
    : params(std::move(params)),
      mat(this->params.nNodes * twistLength, this->params.nNodes * twistLength),
      lhs(this->params.nNodes * twistLength, 1), rhs(this->params.nNodes * twistLength, 1),
      logger(std::cout)
{
    logger.log<LogLevel::info>() << "Initializing solver with " << params.nNodes << " nodes.";
    logger.startTimer("initialization");

    // Initialize data structures
    initCoordinates();
    initSystemMatrix();
    
    logger.endTimer();
    return;
}

StaticSolver::~StaticSolver()
{
}

void StaticSolver::addForce(ScalarType s, Eigen::Ref<const TwistType> force, bool bodyFrame)
{
    // First, transform the force to its nearest node
    CoordType g;
    CoordType fHat = CoordType::Zero();
    ForceListEntry newForce;
    int node = std::round(s / params.ds);
    ScalarType deltaS = s - (node * params.ds);
    se3hat(lhs.getIndexer<twistLength, 1>(BlockIndex::strain, 0).block(lhs.getMat(), node, 0), fHat);
    g = (fHat * deltaS).exp();

    // Calculate the adjoint of the transformation
    SingleMatrixType Adg;
    Adg.setZero();
    Adjoint(g.inverse(), Adg);

    // Fill in the entry and move it over
    newForce.force = Adg.transpose() * force;
    newForce.node = node;
    newForce.bodyFrame = bodyFrame;
    appliedForces.emplace_back(std::move(newForce));

    return;
}

void StaticSolver::addPointForce(ScalarType s, Eigen::Ref<const TwistType> force, bool bodyFrame)
{
    // If the force is being placed at the end, we need to double it
    // since the end covers only half a segment.
    int node = std::round(s / params.ds);
    addForce(s, force/(params.ds*params.length), bodyFrame);
    
    return;
}

void StaticSolver::addDistributedForce(ScalarType s1, ScalarType s2,
				       Eigen::Ref<const TwistType> force, bool bodyFrame)
{
    // First deal with the leftmost node
    int leftNode = std::round(s1 / params.ds);
    ScalarType sNode = leftNode * params.ds;
    ScalarType width = sNode + params.ds/2 - s1;
    ScalarType mid = s1 + width / 2;
    addForce(mid, force * width / params.ds, bodyFrame);
    //addForce(sNode, force * width / params.ds, bodyFrame);
    
    // Now the rightmost node
    int rightNode = std::round(s2 / params.ds);
    sNode = rightNode * params.ds;
    width = sNode + params.ds/2 - s2;
    mid = s2 - width / 2;
    addForce(mid, force * width / params.ds, bodyFrame);
    
    // Interior nodes
    for (int node = leftNode + 1; node < rightNode; node++)
    {
	addForce(params.ds * node, force, bodyFrame);
    }
    return;
}

void StaticSolver::clearForces()
{
    appliedForces.clear();
    return;
}

void StaticSolver::addFixedConstraint(int node, Eigen::Ref<const CoordType> g)
{
    // Add the new constraint to the list
    logger.startTimer("addFixedConstraint");
    ConstraintListEntry constraint;
    constraint.active = true;
    constraint.node = node;
    constraint.g = g;

    fixedConstraints.emplace_back(std::move(constraint));

    // The system matrix / vectors need to be updated to reflect the added constraint
    mat.enlargeBlock<BlockDim::row>(BlockIndex::fixedConstraint,
				    mat.blockSize<BlockDim::row>(BlockIndex::fixedConstraint) - 1,
				    twistLength);
    mat.enlargeBlock<BlockDim::col>(BlockIndex::fixedConstraint,
				    mat.blockSize<BlockDim::col>(BlockIndex::fixedConstraint) - 1,
				    twistLength);
    lhs.enlargeBlock<BlockDim::row>(BlockIndex::fixedConstraint,
				    lhs.blockSize<BlockDim::row>(BlockIndex::fixedConstraint) - 1,
				    twistLength);
    rhs.enlargeBlock<BlockDim::row>(BlockIndex::fixedConstraint,
				    rhs.blockSize<BlockDim::row>(BlockIndex::fixedConstraint) - 1,
				    twistLength);
    generateFixedConstraints();
    logger.endTimer();
    return;
}

void StaticSolver::solve()
{
    logger.startTimer("solve");

    // Prepare
    rhs.clearBlock(BlockIndex::twist, 0);
    rhs.clearBlock(BlockIndex::phi, 0);
    rhs.clearBlock(BlockIndex::strain, 0);
    updateStrainAdjoint();
    updateTwistAdjoint();
    updateAppliedForces();
    updateFixedConstraints();

    // Solve
    sparseLuSolver.compute(mat.getMat());
    if (sparseLuSolver.info() != Eigen::ComputationInfo::Success)
    {
	logger.log<LogLevel::fatal>() << "System computation failed!";
    }
    lhs.getMatRef() = sparseLuSolver.solve(rhs.getMat());

    updateCoordinates();
    
    logger.endTimer();
    return;
}

void StaticSolver::initCoordinates()
{
    logger.startTimer("Coordinate initialization");
    // First, add the necessary blocks to the LHS and RHS vectors
    lhs = BlockMatrix<DenseType>(params.nNodes*twistLength, 1);
    lhs.addBlock<BlockDim::row>(mat.blockSize<BlockDim::row>(BlockIndex::twist)); // Phi block
    lhs.addBlock<BlockDim::row>(mat.blockSize<BlockDim::row>(BlockIndex::twist)); // Strain block
    lhs.addBlock<BlockDim::row>(fixedConstraints.size() * twistLength); // Fixed constraint block
    lhs.getMatRef().setZero();
    rhs = BlockMatrix<DenseType>(params.nNodes*twistLength, 1);
    rhs.addBlock<BlockDim::row>(mat.blockSize<BlockDim::row>(BlockIndex::twist));
    rhs.addBlock<BlockDim::row>(mat.blockSize<BlockDim::row>(BlockIndex::twist));
    rhs.addBlock<BlockDim::row>(fixedConstraints.size() * twistLength);
    rhs.getMatRef().setZero();

    // Now we initialize free strain and twists
    fStar = DenseType::Zero(twistLength * params.nNodes, 1);
    gBody = std::vector<CoordType>(params.nNodes);
    Eigen::Vector<ScalarType, 3> p {0, 0, 0};
    Indexer fsInd = Indexer<twistLength, 1>(fStar.rows(), fStar.cols());
    Indexer phiInd = lhs.getIndexer<twistLength, 1>(BlockIndex::phi, 0);
    for (int i = 0; i < params.nNodes; i++)
    {
	// Initial coordinate
	gBody[i] = CoordType::Identity();
	gBody[i].block<3, 1>(0, 3) = p;
	p(0) += params.length / params.nSegments;
    }
    // Compute free strain and phi based on initial geometry
    CoordType fHat;
    for (int i = 0; i < (params.nNodes - 1); i++)
    {
	// Free strain
	fHat = ((gBody[0].inverse()*gBody[i+1]).log() - (gBody[0].inverse()*gBody[i]).log()) / params.ds;
	se3unhat(fHat, fsInd.block(fStar, i, 0));
	// Phi
	se3unhat((gBody[i].inverse()*gBody[i+1]).log(), phiInd.block(lhs.getMatRef(), i+1, 0));
    }
    // Set the free strain of the end to be the same as the previous node
    fsInd.block(fStar, -1, 0) = fsInd.block(fStar, -2, 0);
    // Set initial strains to be equal to the initial free strains
    lhs.copyToBlock(fStar, BlockIndex::strain, 0);

    logger.endTimer();
    return;
}

void StaticSolver::initSystemMatrix()
{
    logger.startTimer("System matrix initialization");
    // Initialize the matrix with the twist dimension
    mat = BlockMatrix<SparseType>(params.nNodes * twistLength, params.nNodes * twistLength);
    // Add the phi dimension
    mat.addBlock<BlockDim::row>(mat.blockSize<BlockDim::row>(0));
    mat.addBlock<BlockDim::col>(mat.blockSize<BlockDim::col>(0));
    // Add the strain dimension
    mat.addBlock<BlockDim::row>(mat.blockSize<BlockDim::row>(0));
    mat.addBlock<BlockDim::col>(mat.blockSize<BlockDim::col>(0));
    // Add the fixed constraint dimension
    mat.addBlock<BlockDim::row>(fixedConstraints.size() * twistLength);
    mat.addBlock<BlockDim::col>(fixedConstraints.size() * twistLength);

    // Initialize the stiffness matrix
    K = SparseType(mat.blockSize<BlockDim::row>(BlockIndex::phi),
		   mat.blockSize<BlockDim::col>(BlockIndex::phi));
    MatrixConstructorList initList;
    initList.add<DiagConstructor>(params.tensileModulus * params.area, 0, 0, K.rows(), twistLength);
    initList.add<DiagConstructor>(params.shearModulus * params.area, 1, 1, K.rows(), twistLength);
    initList.add<DiagConstructor>(params.shearModulus * params.area, 2, 2, K.rows(), twistLength);
    initList.add<DiagConstructor>(params.shearModulus * params.momentX, 3, 3, K.rows(), twistLength);
    initList.add<DiagConstructor>(params.tensileModulus * params.momentY, 4, 4, K.rows(), twistLength);
    initList.add<DiagConstructor>(params.tensileModulus * params.momentZ, 5, 5, K.rows(), twistLength);
    constructMatrix(initList, K);

    // Initialize the derivative matrices
    SparseType Dprime = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.add<DiagConstructor>(-1/params.ds, 0, 0, Dprime.rows(), 1);
    initList.add<DiagConstructor>(1/params.ds, 0, twistLength, Dprime.rows(), 1);
    constructMatrix(initList, Dprime);
    SparseType D = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.add<DiagConstructor>(1, 0, 0, twistLength, 1);
    initList.add<DiagConstructor>(1/(params.ds*params.ds), twistLength, twistLength, D.rows(), 1);
    constructMatrix(initList, D);

    // Initialize the integral matrices
    SparseType Eprime = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.add<DiagConstructor>(0.5, 0, 0, Eprime.rows(), 1);
    initList.add<DiagConstructor>(0.5, 0, twistLength, Eprime.rows(), 1);
    constructMatrix(initList, Eprime);
    
    SparseType E = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.add<DiagConstructor>(1.0/3.0, twistLength, 0, E.rows(), 1);
    initList.add<DiagConstructor>(1.0/6.0, twistLength, twistLength, E.rows(), 1);
    constructMatrix(initList, E);

    // Utility matrices for RHS values
    SparseType T = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.add<DiagConstructor>(1/params.ds, -twistLength, -twistLength, T.rows(), 1);
    constructMatrix(initList, T);

    SparseType R = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.add<DiagConstructor>(1/params.ds, twistLength, 0, R.rows(), 1);
    constructMatrix(initList, R);

    H = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.add<DiagConstructor>(1, 0, 0, H.rows(), 1);
    initList.add<DiagConstructor>(-1, twistLength, 0, H.rows(), 1);
    constructMatrix(initList, H);

    // Begin moving the blocks into the system matrix
    SparseType EprimeInv = SparseType(K.rows(), K.cols());
    SparseType Kinv = SparseType(K.rows(), K.cols());
    I = SparseType(K.rows(), K.cols());
    I.setIdentity();
    sparseLuSolver.analyzePattern(Eprime);
    sparseLuSolver.factorize(Eprime);
    EprimeInv = sparseLuSolver.solve(I);
    sparseLuSolver.analyzePattern(K);
    sparseLuSolver.factorize(K);
    Kinv = sparseLuSolver.solve(I);

    // Theses blocks need to be held onto for future calculations
    EKinv = (E*Kinv).pruned();
    strainEqnFreeStrainPreMultiplier = (K*EprimeInv*T).pruned();

    // These blocks are constant or just have factors added to them, so we
    // can update them in place without needing to save
    SparseType twistTwist = H.pruned();
    SparseType phiPhi = D.pruned();
    SparseType phiStrain = (-R).pruned();
    SparseType strainStrain = (K*EprimeInv*Dprime);

    mat.copyToBlock(twistTwist, BlockIndex::twist, BlockIndex::twist);
    mat.copyToBlock(phiPhi, BlockIndex::phi, BlockIndex::phi);
    mat.copyToBlock(phiStrain, BlockIndex::phi, BlockIndex::strain);
    mat.copyToBlock(strainStrain, BlockIndex::strain, BlockIndex::strain);

    // Initialize the others, but no need to set them yet
    Jksi = SparseType(mat.blockSize<BlockDim::row>(BlockIndex::twist),
		      mat.blockSize<BlockDim::col>(BlockIndex::twist));
    Jphi = SparseType(Jksi);
    Jf = SparseType(Jksi);
    EKinvJf = SparseType(Jksi);
    
    logger.endTimer();
    return;
}

void StaticSolver::updateStrainAdjoint()
{
    // Before regenerating the adjoint matrix, we need to add the current
    // value of AK back into some submatrices
    logger.startTimer("updateStrainAdjoint");
    mat.addToBlock(Jf, BlockIndex::strain, BlockIndex::strain);
    mat.addToBlock(EKinvJf, BlockIndex::phi, BlockIndex::strain);

    // Now we regenerate the strain jacobian
    SingleMatrixType Jtemp;
    SingleMatrixType Kblock = K.block(0, 0, twistLength, twistLength);
    Indexer strainInd = lhs.getIndexer<twistLength, 1>(BlockIndex::strain, 0);
    Eigen::Ref<DenseType> lhsRef = lhs.getMatRef();
    Indexer Jindexer = Indexer<twistLength, twistLength>(Jf.rows(), Jf.cols());
    DenseType initialTerm = DenseType(rhs.blockSize<BlockDim::row>(BlockIndex::strain), 1);
    Indexer initInd = Indexer<twistLength, 1>(initialTerm.rows(), initialTerm.cols());
    TripletList triplets;
    for (int i = 0; i < params.nNodes; i++)
    {
	// Compute the new Jacobian and initial value
	Eigen::Ref<const TwistType> f0 = strainInd.block(lhs.getMat(), i, 0);
	Eigen::Ref<const TwistType> fs0 = initInd.block(fStar, i, 0);
	Eigen::Ref<TwistType> init = initInd.block(initialTerm, i, 0);
	computeStrainJacobian(f0, fs0, Kblock, Jtemp);
	computeStrainInitialTerm(f0, fs0, Kblock, init);
	for (int j = 0; j < twistLength; j++)
	{
	    for (int k = 0; k < twistLength; k++)
	    {
		ScalarType value = Jtemp(j, k);
		if (value != 0)
		{
		    triplets.emplace_back(Jindexer.row(i) + j,
					  Jindexer.col(i) + k, value);
		}
	    }
	}
    }
    Jf.setFromTriplets(triplets.cbegin(), triplets.cend());
    EKinvJf = (EKinv*Jf).pruned();

    // Update blocks that depend on Jf
    mat.subtractFromBlock(Jf, BlockIndex::strain, BlockIndex::strain);
    mat.subtractFromBlock(EKinvJf, BlockIndex::phi, BlockIndex::strain);

    // Update RHS
    DenseType rhsTemp = initialTerm - Jf*lhs.getBlock(BlockIndex::strain, 0);
    rhs.addToBlock(EKinv*rhsTemp, BlockIndex::phi, 0);
    rhsTemp -= strainEqnFreeStrainPreMultiplier*fStar;
    rhs.addToBlock(rhsTemp, BlockIndex::strain, 0);

    logger.endTimer();
    return;
}

void StaticSolver::updateTwistAdjoint()
{
    logger.startTimer("updateTwistAdjoint");
    Indexer phiInd = lhs.getIndexer<twistLength, 1>(BlockIndex::phi, 0);
    Indexer twistInd = lhs.getIndexer<twistLength, 1>(BlockIndex::twist, 0);
    Indexer Jind = Indexer<twistLength, twistLength>(Jphi.rows(), Jphi.cols());
    SingleMatrixType JLksi, JRksi, JRphi;
    TripletList phiTriplets, twistTriplets;
    // Reserve space in the lists to avoid reallocations
    phiTriplets.reserve(params.nNodes*twistLength*twistLength);
    twistTriplets.reserve(2*params.nNodes*twistLength*twistLength);
    // Do the first block
    for (int j = 0; j < twistLength; j++)
    {
        twistTriplets.emplace_back(j, j, 1);
    }
    // Now do the inner blocks
    for (int i = 1; i < params.nNodes; i++)
    {
	Eigen::Ref<const TwistType> ksiPrev = twistInd.block(lhs.getMat(), i-1, 0);
	Eigen::Ref<const TwistType> ksiCurrent = twistInd.block(lhs.getMat(), i, 0);
	Eigen::Ref<const TwistType> phiCurrent = phiInd.block(lhs.getMat(), i, 0);
	// Calculate the individual Jacobian for ksi_(i-1)
	bchJacobian<2, 1>(phiCurrent, ksiCurrent, JLksi);
	bchJacobian<2, 0>(ksiPrev, phiCurrent, JRksi);
	bchJacobian<2, 1>(ksiPrev, phiCurrent, JRphi);
	// Copy the adjoint blocks over to the matrices
	for (int j = 0; j < twistLength; j++)
	{
	    for (int k = 0; k < twistLength; k++)
	    {
		int row = Jind.row(i, j);
		int col = Jind.col(i, k);
		ScalarType value = JLksi(j, k);
		if (value != 0.0)
		{
		    twistTriplets.emplace_back(row, col, value);
		}
		value = JRphi(j, k);
		if (value != 0.0)
		{
		    phiTriplets.emplace_back(row, col, value);
		}
		col -= twistLength;
		value = JRksi(j, k);
		if (value != 0.0)
		{
		    twistTriplets.emplace_back(row, col, -value);
		}
	    }
	}
    }
    // Update system matrix
    Jksi.setFromTriplets(twistTriplets.cbegin(), twistTriplets.cend());
    Jphi.setFromTriplets(phiTriplets.cbegin(), phiTriplets.cend());
    mat.copyToBlock(Jksi, BlockIndex::twist, BlockIndex::twist);
    mat.copyToBlock(-Jphi, BlockIndex::twist, BlockIndex::phi);
    // Update RHS
    rhs.copyToBlock(Jksi*lhs.getBlock(BlockIndex::twist, 0) - Jphi*lhs.getBlock(BlockIndex::phi, 0),
		    BlockIndex::twist, 0);
    
    logger.endTimer();
    return;
}

void StaticSolver::updateAppliedForces()
{
    // First, generate a vector of the applied forces
    DenseType forces = DenseType::Zero(rhs.blockSize<BlockDim::row>(BlockIndex::strain), 1);
    Indexer indexer = Indexer<twistLength, 1>(forces.rows(), forces.cols());

    // Iterate through the force list and add all of the forces
    CoordType g;
    SingleMatrixType Adg;
    for (auto it = appliedForces.cbegin(); it != appliedForces.cend(); ++it)
    {	
	if (it->bodyFrame) // Force can be applied directly with no transformation
	{
	    indexer.block(forces, it->node, 0) += it->force;
	}
	else // Need to rotate from the spatial frame into the body frame
	{
	    // Construct a coordinate transformation which only rotates
	    g.setIdentity();
	    g.block<3, 3>(0, 0) = gBody[it->node].block<3, 3>(0, 0);
	    Adg.setZero();
	    Adjoint(g, Adg);
	    indexer.block(forces, it->node, 0) += Adg.transpose() * it->force;
	}
    }
    // Now, the forces are added to the RHS
    rhs.subtractFromBlock(EKinv*forces, BlockIndex::phi, 0);
    rhs.subtractFromBlock(forces, BlockIndex::strain, 0);
    
    return;
}

void StaticSolver::updateFixedConstraints()
{
    // For each active constraint, we simply need to enforce that the twist
    // is equal to the current discrepancy from the actual coordinate
    Indexer rhsInd = rhs.getIndexer<twistLength, 1>(BlockIndex::fixedConstraint, 0);
    CoordType g;
    for (int i = 0; i < fixedConstraints.size(); i++)
    {
	if (fixedConstraints[i].active)
	{
	    se3inv(gBody[fixedConstraints[i].node], g);
	    logger.log() << "Constraint " << i << " desired location:\n" << fixedConstraints[i].g;
	    
	    logger.log() << "Constraint " << i << " current offset:\n" << g*fixedConstraints[i].g;
	    g = (g*fixedConstraints[i].g).log();
	    se3unhat(g, rhsInd.block(rhs.getMatRef(), i, 0));
	    logger.log() << "New twist:\n" << rhsInd.block(rhs.getMat(), i, 0);
	}
	else
	{
	    rhsInd.block(rhs.getMatRef(), i, 0).setZero(); 
	}
    }
    return;
}

void StaticSolver::generateFixedConstraints()
{
    // We start by clearing the row and column of the block matrix associated with
    // the fixed constraints.
    mat.clearDim<BlockDim::row>(BlockIndex::fixedConstraint);
    mat.clearDim<BlockDim::col>(BlockIndex::fixedConstraint);

    // Generate the fixed constraint location matrix and deactivation matrix
    Lf = SparseType(mat.blockSize<BlockDim::row>(BlockIndex::phi),
			       mat.blockSize<BlockDim::col>(BlockIndex::fixedConstraint));
    SparseType sigmaF = SparseType(mat.blockSize<BlockDim::row>(BlockIndex::fixedConstraint),
				   mat.blockSize<BlockDim::col>(BlockIndex::fixedConstraint));
    Indexer ind = Indexer<twistLength, twistLength>(Lf.rows(), Lf.cols());
    for (int i = 0; i < fixedConstraints.size(); i++)
    {
	if (fixedConstraints[i].active)
	{
	    // Set up the matrix block
	    for (int j = 0; j < twistLength; j++)
	    {
		Lf.coeffRef(ind.row(fixedConstraints[i].node, j),
			    ind.col(i, j)) = 1;
	    }
	}
	else
	{
	    for (int j = 0; j < twistLength; j++)
	    {
		sigmaF.coeffRef(ind.row(i, j), ind.col(i, j));
	    }
	}
    }
    // Add the values into the system matrix
    mat.copyToBlock(EKinv*Lf, BlockIndex::phi, BlockIndex::fixedConstraint);
    mat.copyToBlock(Lf.transpose(), BlockIndex::fixedConstraint, BlockIndex::twist);
    mat.copyToBlock(sigmaF, BlockIndex::fixedConstraint, BlockIndex::fixedConstraint);
    mat.copyToBlock(Lf, BlockIndex::strain, BlockIndex::fixedConstraint);

    return;
}

void StaticSolver::updateCoordinates()
{
    CoordType gInv, expPhiHat, ksiHat;
    Indexer phiInd = lhs.getIndexer<twistLength, 1>(BlockIndex::phi, 0);
    Indexer twistInd = lhs.getIndexer<twistLength, 1>(BlockIndex::twist, 0);
    for (int i = 1; i < params.nNodes; i++)
    {
	// Save the inverse of the coordinate from the previous iteration
	se3inv(gBody[i], gInv);
	// Update the coordinate
	se3hat(phiInd.block(lhs.getMat(), i, 0), expPhiHat);
	expPhiHat = expPhiHat.exp();
	gBody[i] = gBody[i-1]*expPhiHat;
	// Update ksi
	//se3unhat((gInv*gBody[i]).log(), twistInd.block(lhs.getMatRef(), i, 0));
	twistInd.block(lhs.getMatRef(), i, 0).setZero();
    }
    return;
}

} // namespace cossolve
