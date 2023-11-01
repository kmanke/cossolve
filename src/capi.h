/* This file defines a C API for cossolve.
 *
 */

#ifndef COSSOLVE_CAPI_H
#define COSSOLVE_CAPI_H

#include "Solver.h"

extern "C" {

    // Hide these in an anonymous namespace to avoid polluting the global namespace
    typedef cossolve::Solver::ScalarType ScalarType;
    typedef void* SolverHandle;

    // Wrapper to construct a Solver
    SolverHandle cossolve_createSolver(int nNodes, ScalarType tensileModulus, ScalarType poissonRatio,
				       ScalarType shearModulus, ScalarType momentX, ScalarType momentY,
				       ScalarType momentZ, ScalarType area, ScalarType length,
				       ScalarType linearDensity);
    // Wrapper to destroy a Solver
    void cossolve_deleteSolver(SolverHandle handle);
    // Fills the passed ScalarType array with the strain vector from this solver
    void cossolve_getStrains(SolverHandle handle, ScalarType* outArray);
    void cossolve_getCoords(SolverHandle handle, ScalarType* outArray);
    int cossolve_getNodeCount(SolverHandle handle);

    // Wrapper for Solver::addForce
    void cossolve_Solver_addForce(SolverHandle handle, ScalarType s, ScalarType* force, bool bodyFrame);

    // Wrapper for Solver::solveStrains
    void cossolve_Solver_solveStrains(SolverHandle handle);

    // Wrapper for Solver::solveCoords
    void cossolve_Solver_solveCoords(SolverHandle handle);
}

#endif // COSSOLVE_CAPI_H
