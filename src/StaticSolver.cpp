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
    if (node == 0 || node == (params.nNodes - 1))
    {
	addForce(s, 2*force/(params.ds*params.length));
    }
    else
    {
	addForce(s, force/(params.ds*params.length));
    }
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
    //addForce(sNode, force * width / params.ds, bodyFrame);
    
    // Interior nodes
    for (int node = leftNode + 1; node < rightNode; node++)
    {
	addForce(params.ds * node, force, bodyFrame);
    }
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
    updateFixedConstraints();
    logger.endTimer();
    return;
}

void StaticSolver::solve()
{
    logger.startTimer("solve");
    rhs.clearBlock(BlockIndex::twist, 0);
    rhs.clearBlock(BlockIndex::strain, 0);
    updateAdjoint();
    updateAppliedForces();
    
    sparseLuSolver.compute(mat.getMat());
    lhs.getMatRef() = sparseLuSolver.solve(rhs.getMat());

    // Update the coordinates
/*    Indexer ind = lhs.getIndexer<twistLength, 1>(BlockIndex::twist, 0);
    for (int i = 0; i < params.nNodes; i++)
    {
	gBody[i].setZero();
	se3hat(ind.block(lhs.getMat(), i, 0), gBody[i]);
	gBody[i] = gBody[i].exp();
	}*/
    Indexer ind = lhs.getIndexer<twistLength, 1>(BlockIndex::strain, 0);
    for (int i = 1; i < params.nNodes; i++)
    {
	CoordType fHat = CoordType::Zero();
	se3hat(ind.block(lhs.getMat(), i-1, 0), fHat);
	gBody[i] = gBody[i-1]*(fHat * params.ds).exp();
    }
    logger.endTimer();
    return;
}

void StaticSolver::initCoordinates()
{
    logger.startTimer("Coordinate initialization");
    // First, add the necessary blocks to the LHS and RHS vectors
    lhs.addBlock<BlockDim::row>(mat.blockSize<BlockDim::row>(BlockIndex::twist));
    lhs.addBlock<BlockDim::row>(fixedConstraints.size() * twistLength);
    rhs.addBlock<BlockDim::row>(mat.blockSize<BlockDim::row>(BlockIndex::twist));
    rhs.addBlock<BlockDim::row>(fixedConstraints.size() * twistLength);

    // Now we initialize free strain and coordinate transformations
    fStar = DenseType::Zero(twistLength * params.nNodes, 1);
    Indexer indexer = Indexer<twistLength, 1>(fStar.rows(), fStar.cols());
    gBody = std::vector<CoordType>(params.nNodes);
    Eigen::Vector<ScalarType, 3> p {0, 0, 0};
    for (auto gBodyIt = gBody.begin(); gBodyIt != gBody.end(); ++gBodyIt)
    {
	*gBodyIt = CoordType::Identity();
	(*gBodyIt).block<3, 1>(0, 3) = p;

	p(0) += params.length / params.nSegments;
    }
    // Compute free strains based on initial geometry
    CoordType fHat;
    for (int i = 0; i < (params.nNodes - 1); i++)
    {
	fHat = ((gBody[0].inverse()*gBody[i+1]).log() - (gBody[0].inverse()*gBody[i]).log()) / params.ds;
	se3unhat(fHat, indexer.block(fStar, i, 0));
    }
    // Set the free strain of the end to be the same as the previous node
    indexer.block(fStar, -1, 0) = indexer.block(fStar, -2, 0);
    // Set initial strains to be equal to the initial free strains
    lhs.copyToBlock(fStar, BlockIndex::strain, 0);
    logger.endTimer();
    return;
}

void StaticSolver::initSystemMatrix()
{
    logger.startTimer("System matrix initialization");
    // Initialize the matrix and add blocks
    mat = BlockMatrix<SparseType>(params.nNodes * twistLength, params.nNodes * twistLength);
    // Add the strain dimension
    mat.addBlock<BlockDim::row>(mat.blockSize<BlockDim::row>(0));
    mat.addBlock<BlockDim::col>(mat.blockSize<BlockDim::col>(0));
    // Add the fixed constraint dimension
    mat.addBlock<BlockDim::row>(fixedConstraints.size() * twistLength);
    mat.addBlock<BlockDim::col>(fixedConstraints.size() * twistLength);

    // Initialize the stiffness matrix
    K = SparseType(mat.blockSize<BlockDim::row>(BlockIndex::twist),
			      mat.blockSize<BlockDim::col>(BlockIndex::twist));
    MatrixConstructorList initList;
    initList.emplace_back(std::make_unique<DiagConstructor>(params.tensileModulus * params.area, 0, 0, K.rows(), twistLength));
    initList.emplace_back(std::make_unique<DiagConstructor>(params.shearModulus * params.area, 1, 1, K.rows(), twistLength));
    initList.emplace_back(std::make_unique<DiagConstructor>(params.shearModulus * params.area, 2, 2, K.rows(), twistLength));
    initList.emplace_back(std::make_unique<DiagConstructor>(params.shearModulus * params.momentX, 3, 3, K.rows(), twistLength));
    initList.emplace_back(std::make_unique<DiagConstructor>(params.tensileModulus * params.momentY, 4, 4, K.rows(), twistLength));
    initList.emplace_back(std::make_unique<DiagConstructor>(params.tensileModulus * params.momentZ, 5, 5, K.rows(), twistLength));
    constructMatrix(initList, K);

    // Initialize the derivative matrices
    SparseType Dprime = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.emplace_back(std::make_unique<DiagConstructor>(1, 0, 0, twistLength, 1));
    initList.emplace_back(std::make_unique<DiagConstructor>(-1/params.ds, twistLength, 0, K.rows(), 1));
    initList.emplace_back(std::make_unique<DiagConstructor>(1/params.ds, twistLength, twistLength, K.rows(), 1));
    constructMatrix(initList, Dprime);
    SparseType D = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.emplace_back(std::make_unique<DiagConstructor>(-1/(params.ds*params.ds), twistLength, 0, K.rows(), 1));
    initList.emplace_back(std::make_unique<DiagConstructor>(1/(params.ds*params.ds), twistLength, twistLength, K.rows(), 1));
    constructMatrix(initList, D);

    // Initialize the integral matrices
    SparseType Eprime = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.emplace_back(std::make_unique<DiagConstructor>(1, 0, 0, twistLength, 1));
    initList.emplace_back(std::make_unique<DiagConstructor>(0.5, twistLength, 0, Eprime.rows(), 1));
    initList.emplace_back(std::make_unique<DiagConstructor>(0.5, twistLength, twistLength, Eprime.rows(), 1));
    constructMatrix(initList, Eprime);
    SparseType E = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.emplace_back(std::make_unique<RowConstructor>(2, 0, twistLength, -twistLength, twistLength));
    initList.emplace_back(std::make_unique<RowConstructor>(2, 1, twistLength + 1, -twistLength, twistLength));
    initList.emplace_back(std::make_unique<RowConstructor>(2, 2, twistLength + 2, -twistLength, twistLength));
    initList.emplace_back(std::make_unique<RowConstructor>(2, 3, twistLength + 3, -twistLength, twistLength));
    initList.emplace_back(std::make_unique<RowConstructor>(2, 4, twistLength + 4, -twistLength, twistLength));
    initList.emplace_back(std::make_unique<RowConstructor>(2, 5, twistLength + 5, -twistLength, twistLength));
    initList.emplace_back(std::make_unique<DiagConstructor>(1, 0, 0, twistLength, 1));
    initList.emplace_back(std::make_unique<DiagConstructor>(1, 0, -twistLength, twistLength, 1));
    initList.emplace_back(std::make_unique<DiagConstructor>(1.0/3.0, twistLength, 0, Eprime.rows(), 1));
    initList.emplace_back(std::make_unique<DiagConstructor>(1.0/6.0, twistLength, twistLength, Eprime.rows(), 1));    
    constructMatrix(initList, E);

    // Utility matrices for RHS values
    SparseType T = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.emplace_back(std::make_unique<DiagConstructor>(1, 0, 0, twistLength, 1));
    constructMatrix(initList, T);

    SparseType R = SparseType(K.rows(), K.cols());
    initList.clear();
    initList.emplace_back(std::make_unique<DiagConstructor>(1/params.ds, twistLength, 0, R.rows(), 1));
    constructMatrix(initList, R);

    // Begin moving the blocks into the system matrix
    SparseType EprimeInv = SparseType(K.rows(), K.cols());
    SparseType Kinv = SparseType(K.rows(), K.cols());
    SparseType Einv = SparseType(K.rows(), K.cols());
    SparseType I = SparseType(K.rows(), K.cols());
    I.setIdentity();
    sparseLuSolver.compute(Eprime);
    EprimeInv = sparseLuSolver.solve(I);
    sparseLuSolver.compute(K);
    Kinv = sparseLuSolver.solve(I);
    sparseLuSolver.compute(E);
    Einv = sparseLuSolver.solve(I);

    // Theses blocks need to be held onto for future calculations
    strainEqnForcePreMultiplier = (-K*EprimeInv*T*Kinv + I).pruned();
    strainEqnFreeStrainPreMultiplier = (K*EprimeInv*T).pruned();

    // These blocks are constant or just have AK added to them, so we
    // can update them in place without needing to save
    SparseType twistTwist = (K*Einv*D).pruned();
    SparseType twistStrain = (-K*Einv*R).pruned();
    SparseType strainStrain = (K*EprimeInv*Dprime);
    
    mat.copyToBlock(twistTwist, BlockIndex::twist, BlockIndex::twist);
    mat.copyToBlock(twistStrain, BlockIndex::twist, BlockIndex::strain);
    mat.copyToBlock(strainStrain, BlockIndex::strain, BlockIndex::strain);

//    logger.log() << '\n' << DenseType(D);
    //  logger.log() << '\n' << DenseType(E);
    //logger.log() << '\n' << DenseType(R);
    
    logger.endTimer();
    return;
}

void StaticSolver::updateAdjoint()
{
    // Before regenerating the adjoint matrix, we need to add the current
    // value of AK back into some submatrices
    mat.addToBlock(AK, BlockIndex::strain, BlockIndex::strain);
    mat.addToBlock(AK, BlockIndex::twist, BlockIndex::strain);

    // Now we regenerate the adjoint matrix
    SingleMatrixType adf;
    Indexer strainIndexer = lhs.getIndexer<twistLength, 1>(BlockIndex::strain, 0);
    Eigen::Ref<DenseType> lhsRef = lhs.getMatRef();
    SparseType A = SparseType(mat.blockSize<BlockDim::row>(BlockIndex::strain),
			      mat.blockSize<BlockDim::col>(BlockIndex::strain));
    Indexer Aindexer = Indexer<twistLength, twistLength>(A.rows(), A.cols());
    TripletList triplets;
    for (int i = 0; i < params.nNodes; i++)
    {
	adf.setZero();
	adjoint(strainIndexer.block(lhsRef, i, 0), adf);
	adf.transposeInPlace();
	for (int j = 0; j < twistLength; j++)
	{
	    for (int k = 0; k < twistLength; k++)
	    {
		if (adf(j, k) != 0)
		{
		    triplets.emplace_back(Aindexer.row(i) + j,
					  Aindexer.col(i) + k, adf(j, k));
		}
	    }
	}
    }
    A.setFromTriplets(triplets.cbegin(), triplets.cend());
    AK = (A * K).pruned();

    // Update blocks that depend on AK
    mat.subtractFromBlock(AK, BlockIndex::strain, BlockIndex::strain);
    mat.subtractFromBlock(AK, BlockIndex::twist, BlockIndex::strain);

    // Update RHS
    DenseType rhsTemp = AK*fStar;
    rhs.subtractFromBlock(rhsTemp, BlockIndex::twist, 0);
    rhsTemp -= strainEqnFreeStrainPreMultiplier * fStar;
    rhs.subtractFromBlock(rhsTemp, BlockIndex::strain, 0);
    
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
    rhs.subtractFromBlock(forces, BlockIndex::twist, 0);
    forces = strainEqnForcePreMultiplier * forces;
    rhs.subtractFromBlock(forces, BlockIndex::strain, 0);

    return;
}

void StaticSolver::updateFixedConstraints()
{
    // We start by clearing the row and column of the block matrix associated with
    // the fixed constraints.
    mat.clearDim<BlockDim::row>(BlockIndex::fixedConstraint);
    mat.clearDim<BlockDim::col>(BlockIndex::fixedConstraint);

    // Generate the fixed constraint location matrix and deactivation matrix
    SparseType Lf = SparseType(mat.blockSize<BlockDim::row>(BlockIndex::twist),
			       mat.blockSize<BlockDim::col>(BlockIndex::fixedConstraint));
    SparseType sigmaF = SparseType(mat.blockSize<BlockDim::row>(BlockIndex::fixedConstraint),
				   mat.blockSize<BlockDim::col>(BlockIndex::fixedConstraint));
    Indexer ind = Indexer<twistLength, twistLength>(Lf.rows(), Lf.cols());
    Indexer rhsInd = rhs.getIndexer<twistLength, 1>(BlockIndex::fixedConstraint, 0);
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
	    // Set up the RHS
	    se3unhat(fixedConstraints[i].g.log(), rhsInd.block(rhs.getMatRef(), i, 0));
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
    mat.copyToBlock(Lf, BlockIndex::twist, BlockIndex::fixedConstraint);
    mat.copyToBlock(Lf.transpose(), BlockIndex::fixedConstraint, BlockIndex::twist);
    mat.copyToBlock(sigmaF, BlockIndex::fixedConstraint, BlockIndex::fixedConstraint);
    Lf = strainEqnForcePreMultiplier * Lf;
    mat.copyToBlock(Lf, BlockIndex::strain, BlockIndex::fixedConstraint);

    return;
}

} // namespace cossolve
