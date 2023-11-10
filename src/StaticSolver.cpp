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

#include <eigen3/unsupported/Eigen/MatrixFunctions>

#include <iostream>

namespace cossolve {

StaticSolver::StaticSolver(SolverParameters&& params)
    : params(std::move(params)), sysMat(this->params, fixedConstraints),
      sysVec(this->params, appliedForces, fixedConstraints),
      logger(std::cout)
{
    logger.log<LogLevel::info>() << "Initializing solver with " << params.nNodes() << " nodes.";
    logger.startTimer("initialization");

    // Initialize data structures
    initCoordinates();
    TwistType stiffnessDiag
	{ params.tensileModulus() * params.area(), params.shearModulus() * params.area(),
	  params.shearModulus() * params.area(), params.shearModulus() * params.momentX(),
	  params.tensileModulus() * params.momentY(), params.tensileModulus() * params.momentZ() };
    DerivativeMatrixCoefficients derivCoeffs
	{ -1/(2*params.ds()*params.ds()), 1/(2*params.ds()*params.ds()),    // First node coefficients
	  0, -1/(2*params.ds()*params.ds()), 1/(2*params.ds()*params.ds()), // Inner node coefficients
	  0, 0 };                       // Last node coefficients
    IntegralMatrixCoefficients intCoeffs
	{ 2.0 / 3, 1.0/3,        // First node coefficients
	  1, 2, 2, 5.0/3, 1.0/3, // Inner node coefficients
	  1, 2, 2, 1 };          // Last node coefficients
    sysMat.initStrainMatrix(stiffnessDiag, derivCoeffs, intCoeffs, sparseLuSolver);
    sysVec.initAppliedForceConstants(sysMat.getEinv(), fStar);
    sysVec.initStrains(fStar);

    CoordType(g) = CoordType::Identity();
    g(0, 3) = 0.25;
    addFixedConstraint(1, g);

    logger.log() << '\n' << Eigen::MatrixXd(sysMat.getMat());
    logger.log() << '\n' << sysVec.getLhs();
    logger.log() << '\n' << sysVec.getRhs();
    
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
    ForceListEntry newForce;
    g.setZero();
    int node = intermediateTransform(s, g);

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
    int node = std::round(s / params.ds());
    if (node == 0 || node == (params.nNodes() - 1))
    {
	addForce(s, 2*force/(params.ds()*params.length()));
    }
    else
    {
	addForce(s, force/(params.ds()*params.length()));
    }
    return;
}

void StaticSolver::addDistributedForce(ScalarType s1, ScalarType s2,
				       Eigen::Ref<const TwistType> force, bool bodyFrame)
{
    // First deal with the leftmost node
    int leftNode = std::round(s1 / params.ds());
    ScalarType sNode = leftNode * params.ds();
    ScalarType width = sNode + params.ds()/2 - s1;
    ScalarType mid = s1 + width / 2;
    addForce(mid, force * width / params.ds(), bodyFrame);
    //addForce(sNode, force * width / params.ds(), bodyFrame);
    
    // Now the rightmost node
    int rightNode = std::round(s2 / params.ds());
    sNode = rightNode * params.ds();
    width = sNode + params.ds()/2 - s2;
    mid = s2 - width / 2;
    addForce(mid, force * width / params.ds(), bodyFrame);
    //addForce(sNode, force * width / params.ds(), bodyFrame);
    
    // Interior nodes
    for (int node = leftNode + 1; node < rightNode; node++)
    {
	addForce(params.ds() * node, force, bodyFrame);
    }

    return;
}

void StaticSolver::addFixedConstraint(int node, Eigen::Ref<const CoordType> g)
{
    ConstraintListEntry constraint;
    constraint.active = true;
    constraint.node = node;
    constraint.g = g;

    fixedConstraints.emplace_back(std::move(constraint));

    // System needs to be updated with the new constraints
    sysVec.addFixedConstraint(fixedConstraints.size() - 1);
    sysMat.addFixedConstraint(fixedConstraints.size() - 1);
    sysVec.updateFixedConstraints();
    sysMat.updateFixedConstraints();
    
    return;
}

void StaticSolver::timeStep()
{
//    applyForces(1);
//    applyContactForces();
//    computeStrains();
    //   computeDeformation();

    return;
}

void StaticSolver::initCoordinates()
{
    fStar = VectorType::Zero(twistLength * params.nNodes());
    gBody = std::vector<CoordType>(params.nNodes());

    // Initialize coordinate transformations
    Eigen::Vector<ScalarType, 3> p {0, 0, 0};
    for (auto gBodyIt = gBody.begin(); gBodyIt != gBody.end(); ++gBodyIt)
    {
	*gBodyIt = CoordType::Identity();
	(*gBodyIt).block<3, 1>(0, 3) = p;

	p(0) += params.length() / params.nSegments();
    }

    // Compute free strains based on initial geometry
    CoordType fHat;
    for (int i = 0; i < (params.nNodes() - 1); i++)
    {
	fHat = ((gBody[0].inverse()*gBody[i+1]).log() - (gBody[0].inverse()*gBody[i]).log()) / params.ds();
	CoordType dg = (gBody[i+1] - gBody[i]) / params.ds();
	//se3unhat(gBody[i].inverse()*dg, vectorRef(fStar, i, twistLength));
	se3unhat(fHat, vectorRef(fStar, i, twistLength));
    }
    // Set the free strain of the end to be the same as the previous node
    vectorRef(fStar, -1, twistLength) = vectorRef(fStar, -2, twistLength);
    return;
}

void StaticSolver::generateAdjointMatrix()
{
    SingleMatrixType adf;
    for (int i = 0; i < params.nNodes(); i++)
    {
	adf.setZero();
//	adjoint(nodeVectorConst(strains, i), adf);
	adf.transposeInPlace();
	for (int j = 0; j < Sizes::entriesPerVector; j++)
	{
	    for (int k = 0; k < Sizes::entriesPerVector; k++)
	    {
//		Af.coeffRef(nodeIndex(i)+j, nodeIndex(i)+k) = adf(j, k);
	    }
	}
    }
    return;
}

int StaticSolver::intermediateTransform(ScalarType s, Eigen::Ref<CoordType> g) const
{
    // First, find the nearest node to s
    int nearestNode = std::round(s / params.ds());
    ScalarType sk = nearestNode * params.ds();
    
    // Compute the strain integral
    TwistType integratedStrain = TwistType::Zero();
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

void StaticSolver::intermediateStrainIntegral(ScalarType s1, ScalarType s2, Eigen::Ref<TwistType> result) const
{
    // Calculate the first node before s1 and the first node before s2
    int node1 = std::floor(s1 / params.ds());
    int node2 = std::floor(s2 / params.ds());

    // First, integrate from node1 to node2
    nodalStrainIntegral(node1, node2, result);

    // Now we need to fix up the start and the end
    relativeStrainIntegral(node1, s1 - node1*params.ds(), result, true); // Subtract the first part
    relativeStrainIntegral(node2, s2 - node2*params.ds(), result, false); // Add the end part

    return;
}

void StaticSolver::nodalStrainIntegral(int node1, int node2, Eigen::Ref<TwistType> result) const
{
    VectorType strains(twistLength * (params.nNodes() + 1));
/*    while (node1 < node2)
    {
	if (node1 == 0) // First node
	{
	    Eigen::Ref<const TwistType> f0 = nodeVectorConst(strains, 0);//strains.block<6, 1>(0, 0);
	    Eigen::Ref<const TwistType> f1 = nodeVectorConst(strains, 1);
	    Eigen::Ref<const TwistType> f2 = nodeVectorConst(strains, 2);
	    result += params.ds()*(5.0/12.0*f0 + 8.0/12.0*f1 - 1.0/12.0*f2);
	}
	else // General case
	{
	    Eigen::Ref<const TwistType> fprev = nodeVectorConst(strains, node1 - 1);
	    Eigen::Ref<const TwistType> fk = nodeVectorConst(strains, node1);
	    Eigen::Ref<const TwistType> fnext = nodeVectorConst(strains, node1 + 1);
	    Eigen::Ref<const TwistType> fnext2 = nodeVectorConst(strains, node1 + 2);
	    result += params.ds()*(1.0/6.0*(fnext - fprev) + 11.0/12.0*fk + 1.0/12.0*fnext2);
	}
	node1++;
	}*/
    return;
}

void StaticSolver::relativeStrainIntegral(int node, ScalarType x, Eigen::Ref<TwistType> result, bool subtract) const
{
    /*
    VectorType strains(twistLength * (params.nNodes() + 1));
    if (node == 0) // Special case for the first segment
    {
	Eigen::Ref<const TwistType> f0 = nodeVectorConst(strains, 0);
	Eigen::Ref<const TwistType> f1 = nodeVectorConst(strains, 1);
	Eigen::Ref<const TwistType> f2 = nodeVectorConst(strains, 2);

	result += (subtract ? -1 : 1) * (x*f0 + x*x/(2*params.ds())*(f1 - f0) + x*x*x/(12*params.ds()*params.ds())*(f2 + 2*f1 - 3*f0));
    }
    else if (node < (params.nNodes() - 1)) // General case for middle segments
    {
	Eigen::Ref<const TwistType> fprev = nodeVectorConst(strains, node - 1);
	Eigen::Ref<const TwistType> fk = nodeVectorConst(strains, node);
	Eigen::Ref<const TwistType> fnext = nodeVectorConst(strains, node + 1);
	Eigen::Ref<const TwistType> fnext2 = nodeVectorConst(strains, node + 2);
	result += (subtract ? -1 : 1) * (x*fk + x*x/(4*params.ds())*(fnext - fprev) + x*x*x/(12*params.ds()*params.ds())*(fnext2 - fk + fnext - fprev));
    }
    // If node >= nNodes - 1, then we're at the tip, so do nothing at all */
    return;
}

void StaticSolver::computeStrains()
{
    logger.startTimer("Solve Strain");
    // Reconstruct the adjoint matrices
    logger.startTimer("Adjoint Matrix");
    generateAdjointMatrix();
//    AfK = (Af * K).pruned();
    logger.endTimer();

    // The last vector in our forces vector holds the boundary condition (strain at the imaginary (n+1)th node)
    /*  SingleMatrixType invK = SingleMatrixType(K.topLeftCorner(Sizes::entriesPerVector,
							     Sizes::entriesPerVector)).inverse();
    nodeVector(forces, -1) = -(nodeVector(strains, -2) - ds/2*invK*(nodeVector(forces, -2)));
    */    
/*    // Compute
    logger.startTimer("compute");
    sparseLuSolver.compute((EinvKD - AfK).pruned());
    logger.endTimer();
    logger.startTimer("compute RHS");
    VectorType rhs = -(AfK * freeStrains) - forces;
    logger.endTimer();
    strains = sparseLuSolver.solve(rhs);
    logger.endTimer();*/
    
    return;
}

void StaticSolver::computeDeformation()
{
    TwistType integratedStrain;
    CoordType integratedStrainHat = CoordType::Zero();
/*    for (int node = baseNode - 1; node >= 0; node--)
    {
	integratedStrain.setZero();
	nodalStrainIntegral(node, node + 1, integratedStrain);
	se3hat(integratedStrain, integratedStrainHat);
	gBody[node] = gBody[node+1]*(-integratedStrainHat.exp());
    }
    for (int node = baseNode + 1; node < params.nNodes(); node++)
    {
	integratedStrain.setZero();
	nodalStrainIntegral(node - 1, node, integratedStrain);
	se3hat(integratedStrain, integratedStrainHat);
	gBody[node] = gBody[node-1]*integratedStrainHat.exp();
	}*/
    return;
}

} // namespace cossolve
