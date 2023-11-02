/*
 * Implementation of C API functions.
 */

#include "capi.h"
#include "Solver.h"

#include <iostream>

extern "C" {

SolverHandle cossolve_createSolver(int nNodes, ScalarType tensileModulus, ScalarType poissonRatio,
				   ScalarType shearModulus, ScalarType momentX, ScalarType momentY,
				   ScalarType momentZ, ScalarType area, ScalarType length,
				   ScalarType linearDensity)
{
    return new cossolve::Solver(nNodes, tensileModulus, poissonRatio, shearModulus, momentX, momentY,
				momentZ, area, length, linearDensity);
}

void cossolve_deleteSolver(SolverHandle handle)
{
    delete reinterpret_cast<cossolve::Solver*>(handle);
    return;
}

void cossolve_getStrains(SolverHandle handle, ScalarType* outputArray)
{
    auto strains = reinterpret_cast<cossolve::Solver*>(handle)->getStrains();
    memcpy(outputArray, strains.data(), strains.size() * sizeof(ScalarType));
    
    return;
}

void cossolve_getCoords(SolverHandle handle, ScalarType* outputArray)
{
    auto coords = reinterpret_cast<cossolve::Solver*>(handle)->getCoords();
    for (auto it = coords.cbegin(); it != coords.cend(); ++it)
    {
	memcpy(outputArray, (*it).data(), cossolve::Solver::Sizes::byteCountPerCoord);
	outputArray += cossolve::Solver::Sizes::entriesPerCoord;
    }

    return;
}
    
int cossolve_getNodeCount(SolverHandle handle)
{
    return reinterpret_cast<cossolve::Solver*>(handle)->getNodeCount();
}

void cossolve_Solver_addForce(SolverHandle handle, ScalarType s, ScalarType* force, bool bodyFrame)
{
    reinterpret_cast<cossolve::Solver*>(handle)->addForce
	(s, cossolve::Solver::SingleVectorType(force), bodyFrame);

    return;
}

void cossolve_Solver_solveStrains(SolverHandle handle)
{
//    reinterpret_cast<cossolve::Solver*>(handle)->solveStrains();
    reinterpret_cast<cossolve::Solver*>(handle)->timeStep();
    return;
}
    
void cossolve_Solver_solveCoords(SolverHandle handle)
{
    reinterpret_cast<cossolve::Solver*>(handle)->solveCoords();
    return;
}

}
