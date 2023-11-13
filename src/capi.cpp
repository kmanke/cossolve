/*
 * Implementation of C API functions.
 */

#include "capi.h"
#include "StaticSolver.h"
#include "SolverParameters.h"

using namespace cossolve;

extern "C" {

    SolverHandle cossolve_StaticSolver_construct(int nNodes, ScalarType tensileModulus,
						 ScalarType poissonRatio, ScalarType shearModulus,
						 ScalarType momentX, ScalarType momentY, ScalarType momentZ,
						 ScalarType area, ScalarType length, ScalarType linearDensity)
{
    SolverParameters params(tensileModulus, poissonRatio, shearModulus,
			    momentX, momentY, momentZ, area, length, linearDensity, nNodes);
    return new StaticSolver(std::move(params));
}

void cossolve_StaticSolver_delete(SolverHandle handle)
{
    delete reinterpret_cast<StaticSolver*>(handle);
    return;
}

void cossolve_StaticSolver_addPointForce(SolverHandle handle, ScalarType s,
					 ScalarType* force, bool bodyFrame)
{
    reinterpret_cast<StaticSolver*>(handle)->addPointForce(s, TwistType(force), bodyFrame);
    return;
}
    
void cossolve_StaticSolver_addDistributedForce(SolverHandle handle, ScalarType s1, ScalarType s2,
						   ScalarType* force, bool bodyFrame)
{
    reinterpret_cast<StaticSolver*>(handle)->addDistributedForce(s1, s2, TwistType(force), bodyFrame);
    return;
}

void cossolve_StaticSolver_addFixedConstraint(SolverHandle handle, int node, ScalarType* g)
{
    reinterpret_cast<StaticSolver*>(handle)->addFixedConstraint(node, CoordType(g));
    return;
}

void cossolve_StaticSolver_getCoords(SolverHandle handle, ScalarType* outputArray)
{
    const std::vector<CoordType>& coords = reinterpret_cast<StaticSolver*>(handle)->getCoords();
    for (auto it = coords.cbegin(); it != coords.cend(); ++it)
    {
	memcpy(outputArray, it->data(), it->rows() * it->cols() * sizeof(ScalarType));
	outputArray += it->rows() * it->cols();
    }
    return;
}
    
void cossolve_StaticSolver_getStrains(SolverHandle handle, ScalarType* outputArray)
{
    const DenseType& strains = reinterpret_cast<StaticSolver*>(handle)->getStrains();
    memcpy(outputArray, strains.data(), strains.size() * sizeof(ScalarType));
    
    return;
}

    
void cossolve_StaticSolver_getTwists(SolverHandle handle, ScalarType* outputArray)
{
    const DenseType& twists = reinterpret_cast<StaticSolver*>(handle)->getTwists();
    memcpy(outputArray, twists.data(), twists.size() * sizeof(ScalarType));
    
    return;
}
        
void cossolve_StaticSolver_getFixedConstraintForces(SolverHandle handle, ScalarType* outputArray)
{
    const DenseType& forces = reinterpret_cast<StaticSolver*>(handle)->getFixedConstraintForces();
    memcpy(outputArray, forces.data(), forces.size() * sizeof(ScalarType));
    
    return;
}
           
void cossolve_StaticSolver_getSystemMatrix(SolverHandle handle, ScalarType* outputArray)
{
    DenseType mat = DenseType(reinterpret_cast<StaticSolver*>(handle)->getSystemMatrix());
    memcpy(outputArray, mat.data(), mat.rows() * mat.cols() * sizeof(ScalarType));
    
    return;
}

int cossolve_StaticSolver_getSystemRows(SolverHandle handle)
{
    return reinterpret_cast<StaticSolver*>(handle)->getSystemMatrix().rows();
}
    
int cossolve_StaticSolver_getSystemCols(SolverHandle handle)
{
    return reinterpret_cast<StaticSolver*>(handle)->getSystemMatrix().cols();
}
    
void cossolve_StaticSolver_solve(SolverHandle handle)
{
    reinterpret_cast<StaticSolver*>(handle)->solve();
    return;
}

}
