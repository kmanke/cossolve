/* This file defines a C API for cossolve.
 *
 */

#ifndef COSSOLVE_CAPI_H
#define COSSOLVE_CAPI_H

#include "StaticSolver.h"

extern "C" {

    typedef cossolve::ScalarType ScalarType;
    typedef void* SolverHandle;

    SolverHandle cossolve_StaticSolver_construct(
	int nNodes, ScalarType tensileModulus, ScalarType poissonRatio,
	ScalarType shearModulus, ScalarType momentX, ScalarType momentY,
	ScalarType momentZ, ScalarType area, ScalarType length,
	ScalarType linearDensity
	);
    void cossolve_StaticSolver_delete(SolverHandle handle);
    void cossolve_StaticSolver_addPointForce(SolverHandle handle, ScalarType s,
					     ScalarType* force, bool bodyFrame);
    void cossolve_StaticSolver_addDistributedForce(SolverHandle handle,
						   ScalarType s1, ScalarType s2,
						   ScalarType* force, bool bodyFrame);
    void cossolve_StaticSolver_addFixedConstraint(SolverHandle handle, int node, ScalarType* g);
    void cossolve_StaticSolver_getCoords(SolverHandle handle, ScalarType* outArray);
    void cossolve_StaticSolver_getStrains(SolverHandle handle, ScalarType* outArray);
    void cossolve_StaticSolver_getTwists(SolverHandle handle, ScalarType* outArray);
    void cossolve_StaticSolver_getFixedConstraintForces(SolverHandle handle, ScalarType* outArray);
    void cossolve_StaticSolver_getSystemMatrix(SolverHandle handle, ScalarType* outArray);
    int cossolve_StaticSolver_getSystemRows(SolverHandle handle);
    int cossolve_StaticSolver_getSystemCols(SolverHandle handle);
    void cossolve_StaticSolver_solve(SolverHandle handle);
}

#endif // COSSOLVE_CAPI_H
