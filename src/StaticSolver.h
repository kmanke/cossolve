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
#include "SystemMatrix.h"
//#include "SystemVector.h"
#include "SolverParameters.h"
#include "ForceList.h"

#include <eigen3/Eigen/Eigen>

namespace cossolve {

class StaticSolver
{    
public:
    // This struct contains definitions of several commonly used sizes
    struct Sizes
    {
	static constexpr std::size_t entriesPerVector = 6;
	static constexpr std::size_t entriesPerCoord = CoordType::SizeAtCompileTime;
	static constexpr std::size_t entriesPerPoint = 3;
	static constexpr std::size_t byteCountPerSingleVector = entriesPerVector * sizeof(ScalarType);
	static constexpr std::size_t byteCountPerCoord = entriesPerCoord * sizeof(ScalarType);
	static constexpr std::size_t byteCountPerPoint = entriesPerPoint * sizeof(ScalarType);
    };
    
    StaticSolver(SolverParameters&& params);
    ~StaticSolver();

//    inline Eigen::Ref<const VectorType> getStrains() const { return strains.block(0, 0, nodeIndex(-1), 1); }
    inline const std::vector<CoordType>& getCoords() const { return gBody; }
    inline Eigen::Ref<const VectorType> getStrains() const { return VectorType(); }
    inline int getNodeCount() const { return params.nNodes(); }

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
    
    // Solves a single solution time step
    void timeStep();
    
    // Solves the ODE for strains based on the current loads
    void computeStrains();
    
private:
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
    SystemMatrix<DenseType> lhs;
    SystemMatrix<DenseType> rhs;
    VectorType fStar;
    std::vector<CoordType> gBody; // From body to spatial
    
    // System parameter matrices
    SystemMatrix<SparseType> mat; // Full system matrix
    SparseType K; // Stiffness matrix
    SparseType Kinv; // K^-1
    SparseType E; // Twist integral
    SparseType D; // Twist derivative
    SparseType Eprime; // Strain integral
    SparseType Dprime; // Derivative
    SparseType Af; // Adjoint matrix
    SparseType AfK; // Af * K

    // Solvers
    SparseLUSolver<SparseType> sparseLuSolver;
    
    // Initializes the coordinate transformations and free strains to initial values.
    void initCoordinates();

    // Regenerates the strain adjoint matrix.
    void generateAdjointMatrix();
        
    // This function generates a coordinate transformation matrix which transforms
    // from the nearest node to s to s.
    // Input:
    //   s - The point along the rod for which to generate g.
    //   g - The matrix where the output will be stored.
    // Returns:
    //   The index of the nearest node to s.
    int intermediateTransform(ScalarType s, Eigen::Ref<CoordType> g) const;

    // Integrates the strain between any two points on the rod.
    // Input:
    //   s - The point along the rod. Must be in [0, 1].
    //   f - A vector which will be filled with the result.
    void intermediateStrainIntegral(ScalarType s1, ScalarType s2, Eigen::Ref<TwistType> f) const;

    // Integrates the strain between any two nodes on the rod.
    void nodalStrainIntegral(int node1, int node2, Eigen::Ref<TwistType> result) const;

    // Integrates the strain from a node to a point x away from the node (in the positive direction)
    void relativeStrainIntegral(int node, ScalarType x, Eigen::Ref<TwistType> result, bool subtract = false) const;

    // Updates the body coordinate transformations as a result of
    // body deformation.
    void computeDeformation();
};

} // namespace cossolve

#endif // COSSOLVE_STATIC_SOLVER_H
