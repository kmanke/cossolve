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

#ifndef COSSOLVE_SOLVER_H
#define COSSOLVE_SOLVER_H

#include "config.h"
#include "Logger.h"

#include <eigen3/Eigen/Eigen>

namespace cossolve {

class Solver
{
    using ForceEntry = std::tuple<Eigen::Vector<double, 6>, int, bool>;
    using ForceList = std::list<ForceEntry>;
    
public:
    // Type aliases
    using ScalarType = double; // This defines the data type of all numerical quantities used by the solver
    using VectorType = Eigen::VectorX<ScalarType>;
    using SingleVectorType = Eigen::Vector<ScalarType, 6>;
    using MatrixType = Eigen::SparseMatrix<ScalarType>;
    using SingleMatrixType = Eigen::Matrix<ScalarType, 6, 6>;
    using CoordType = Eigen::Matrix<ScalarType, 4, 4>;
    using PointType = Eigen::Vector<ScalarType, 3>;

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
    
    Solver(int nNodes, ScalarType tensileModulus, ScalarType poissonRatio, ScalarType shearModulus,
	   ScalarType momentX, ScalarType momentY, ScalarType momentZ,
	   ScalarType area, ScalarType length, ScalarType linearDensity);
    ~Solver();

    inline Eigen::Ref<const VectorType> getStrains() const { return strains.block(0, 0, nodeIndex(-1), 1); }
    inline const std::vector<CoordType>& getCoords() const { return gBody; }
    inline int getNodeCount() const { return nNodes; }

    // Adds a force to the rod.
    // The force is transformed from its actual position to act
    // on the nearest node.
    // The magnitude of force must be converted to a force per unit length first.
    void addForce(ScalarType s, Eigen::Ref<const SingleVectorType> force, bool bodyFrame = true);

    // Adds a point force to the rod.
    // This is a convenience function which simply scales the input appropriately for
    // the call to addForce.
    void addPointForce(ScalarType s, Eigen::Ref<const SingleVectorType> force, bool bodyFrame = true);
    
    // Adds a distributed force to the rod.
    // The distributed force is transformed into discrete forces
    // distributed approximately a segment length apart.
    void addDistributedForce(ScalarType s1, ScalarType s2,
			     Eigen::Ref<const SingleVectorType> force, bool bodyFrame = true);
    
    // Solves a single solution time step
    void timeStep();
    
    // Solves the ODE for strains based on the current loads
    void computeStrains();
    
private:
    // Useful utilities
    Logger<LogLevel::debug> logger;
    
    // Physical parameters
    const ScalarType tensileModulus;
    const ScalarType poissonRatio;
    const ScalarType shearModulus;
    const ScalarType momentX;
    const ScalarType momentY;
    const ScalarType momentZ;
    const ScalarType area;
    const ScalarType length;
    const ScalarType linearDensity;

    // Applied loads
    // Each list item is a pair with the first item being the force vector
    // and the second a boolean representing whether this force's direction
    // is in the body frame or the spatial frame.
    ForceList externalForces;
    
    // Discretization
    const int nNodes; // The number of nodes to place on the rod
    const int nSegments; // The number of segments along the rod
    const int nDims; // The number of dimensions in the problem.
    const int endIndex; // Index of one-past-the end of a full-length vector
    const ScalarType ds; // The change in arclength parameter per segment
    ScalarType dt; // The time step.
    ScalarType t;
    
    // System state
    VectorType strains;
    VectorType freeStrains;
    VectorType forces;
    CoordType gCentroid; // From centroid to spatial
    std::vector<CoordType> gBody; // From body to spatial
    std::vector<CoordType> gBodyCentroid; // From centroid to body
    std::vector<CoordType> gCentroidBody; // From body to centroid
    SingleVectorType rigidBodyVelocity;
    SingleVectorType rigidBodyAcceleration;
    VectorType nodalVelocities;
    VectorType nodalAccelerations;
    
    // System parameter matrices
    MatrixType K; // Stiffness matrix
    MatrixType M; // Inertia matrix
    MatrixType Af; // Adjoint strain matrix
    MatrixType Av; // Adjoint velocity matrix
    MatrixType D; // Derivative matrix
    MatrixType E; // Integral matrix
    MatrixType AfK; // Af*K (intermediate matrix used in strain calculation)
    MatrixType AvM; // Av*M (intermediate matrix used in strain calculation)
    MatrixType EinvKD; // inv(E)*K*D
    SingleMatrixType Mrigid; // Rigid body inertia tensor

    // Solvers
    SparseLUSolver<MatrixType> sparseLuSolver;
    
    // Helper functions

    // Returns the index of the first element of the specified node.
    // If node is < 0, returns the first element of the specified node
    // counting from the back (starting at the imaginary nth node)
    inline int nodeIndex(int node) const
    {
	return (node >= 0) ? (node * Sizes::entriesPerVector)
	    : ((nDims + node) * Sizes::entriesPerVector);
    }
    
    // Returns a reference to the vector associated with the specified node.
    // If node is < 0, returns the vector counting from the back.
    inline Eigen::Ref<SingleVectorType> nodeVector(Eigen::Ref<VectorType> vector, int node) const
    {
	return vector.block<Sizes::entriesPerVector, 1>(nodeIndex(node), 0);
    }
    inline Eigen::Ref<const SingleVectorType> nodeVectorConst(Eigen::Ref<const VectorType> vector, int node) const
    {
	return vector.block<Sizes::entriesPerVector, 1>(nodeIndex(node), 0);
    }

    // Returns a (const) reference to the point defining the origin
    // of coordinate frame `g`.
    inline Eigen::Ref<PointType> origin(Eigen::Ref<CoordType> g)
    {
	return g.block<Sizes::entriesPerPoint, 1>(0, 3);
    }
    inline Eigen::Ref<const PointType> originConst(Eigen::Ref<const CoordType> g)
    {
	return g.block<Sizes::entriesPerPoint, 1>(0, 3);
    }
    void initStateVectors();

    // The following functions initialize the various system matrices.
    // initSystemMatrices calls each subsequent function to do full initialization.
    void initSystemMatrices();
    void initStiffnessMatrix();
    void initInertiaMatrix();
    void initDerivativeMatrix();
    void initIntegralMatrix();
    void initStrainEqnMatrix();

    // Regenerates the strain adjoint matrix.
    void generateAdjointMatrix();
        
    // This function places forces on the beam.
    // The force is transformed from its actual position s to the nearest node on
    // the rod.
    void applyForces(ScalarType scaleFactor);

    // Calculates contact forces with the ground and applies them to the rod
    void applyContactForces();

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
    void intermediateStrainIntegral(ScalarType s1, ScalarType s2, Eigen::Ref<SingleVectorType> f) const;

    // Integrates the strain between any two nodes on the rod.
    void nodalStrainIntegral(int node1, int node2, Eigen::Ref<SingleVectorType> result) const;

    // Integrates the strain from a node to a point x away from the node (in the positive direction)
    void relativeStrainIntegral(int node, ScalarType x, Eigen::Ref<SingleVectorType> result, bool subtract = false) const;

    // Computes the centroid and mass moment of inertia of the rod based
    // on the current geometry.
    void computeInertia();

    // Computes the rigid body acceleration and velocity of the rod.
    void computeRigidKinematics(ScalarType stepSize);

    // Applies the rigid body velocity to all nodes on the rod.
    void computeNodalVelocities();

    // Applies the rigid body acceleration to all nodes on the rod.
    void computeNodalAccelerations();

    // Updates the body coordinate transformations as a result of
    // rigid body motion.
    void computeTranslation();

    // Updates the body coordinate transformations as a result of
    // body deformation.
    void computeDeformation(int baseNode = 0);
    
    // Based on the current shape, solves for the equilibrium position
    // and forces on the body.
    void solveRigidEquilibrium(ScalarType solverTime);
};

} // namespace cossolve

#endif // COSSOLVE_SOLVER_H
