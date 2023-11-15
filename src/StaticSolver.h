/* This class represents a solver which holds simulation parameters and stores results.
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

#ifndef COSSOLVE_STATIC_SOLVER_H
#define COSSOLVE_STATIC_SOLVER_H

#include "config.h"
#include "CossolveTypes.h"
#include "Logger.h"
#include "BlockMatrix.h"
#include "SolverParameters.h"
#include "ForceList.h"

#include <eigen3/Eigen/Eigen>

namespace cossolve {

class StaticSolver
{    
public:    
    StaticSolver(SolverParameters&& params);
    ~StaticSolver();

    inline const std::vector<CoordType>& getCoords() const { return gBody; }
    inline Eigen::Ref<const VectorType> getStrains() const
	{ return lhs.getBlock(BlockIndex::strain, 0); }
    inline Eigen::Ref<const VectorType> getTwists() const
	{ return lhs.getBlock(BlockIndex::twist, 0); }
    inline Eigen::Ref<const VectorType> getFixedConstraintForces() const
	{ return lhs.getBlock(BlockIndex::fixedConstraint, 0); }
    inline const SparseType& getSystemMatrix() const { return mat.getMat(); }

    // Adds a force to the rod.
    // The force is transformed from its actual position to act
    // on the nearest node.
    // The magnitude of force must be converted to a force per unit length first.
    void addForce(ScalarType s, Eigen::Ref<const TwistType> force, bool bodyFrame = true);

    // Adds a point force to the rod.
    // This is a convenience function which simply scales the input appropriately for
    // the call to addForce.
    void addPointForce(ScalarType s, Eigen::Ref<const TwistType> force, bool bodyFrame = true);
    
    // Adds a distributed force to the rod.
    // The distributed force is transformed into discrete forces
    // distributed approximately a segment length apart.
    void addDistributedForce(ScalarType s1, ScalarType s2,
			     Eigen::Ref<const TwistType> force, bool bodyFrame = true);

    // Adds a fixed constraint to the rod.
    void addFixedConstraint(int node, Eigen::Ref<const CoordType> g);
    
    // Solves the system
    void solve();
    
private:
    struct BlockIndex
    {
	static constexpr int twist = 0;
	static constexpr int phi = 1;
	static constexpr int strain = 2;
	static constexpr int fixedConstraint = 3;
    };
    // Useful utilities
    Logger<LogLevel::debug> logger;

    SolverParameters params;
    
    // Applied loads
    // Each list item is a pair with the first item being the force vector
    // and the second a boolean representing whether this force's direction
    // is in the body frame or the spatial frame.
    ForceList appliedForces;
    ConstraintList fixedConstraints;
    
    // System state
    BlockMatrix<DenseType> lhs;
    BlockMatrix<DenseType> rhs;
    DenseType fStar;
    std::vector<CoordType> gBody; // From body to spatial
    
    // System parameter matrices
    BlockMatrix<SparseType> mat; // Full system matrix

    // These submatrices need to be held onto to update the system matrix
    SparseType I; // 6*nNodes x 6*nNodes identity matrix, we need this somewhat often
    SparseType K; // Stiffness matrix
    SparseType H;
    SparseType Jphi;
    SparseType Jksi;
    SparseType Lf; // Fixed constraint location matrix
    SparseType Jf;
    SparseType EKinv;   // E*K^-1
    SparseType EKinvJf; // E*K^-1*Jf
    SparseType strainEqnFreeStrainPreMultiplier;

    // Solvers
    //Eigen::SparseLU<SparseType> sparseLuSolver;
    SparseLUSolver<SparseType> sparseLuSolver;
    
    // Initializes the coordinate transformations and free strains to initial values.
    void initCoordinates();

    // Initializatin subfunctions
    void initSystemMatrix();
    
    // Regenerates the strain adjoint matrix, then updates all `mat` and `rhs` blocks
    // which depend on it.
    void updateStrainAdjoint();

    // Regenerates the phi and ksi adjoint matrices, then updates all blocks
    // which depend on it.
    void updateTwistAdjoint();

    // Updates the RHS with the current applied forces
    void updateAppliedForces();

    // Updates the system matrix to reflect changed fixed constraints.
    // This only needs to be called when fixed constraints change. There is no need
    // to call it every iteration.
    void updateFixedConstraints();
};

} // namespace cossolve

#endif // COSSOLVE_STATIC_SOLVER_H
