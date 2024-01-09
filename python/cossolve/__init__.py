# Python interface for cossolve.
# Copyright (C) 2023 Kyle Manke
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import ctypes

# Import the shared library
lib_path = '/home/kmanke/dev/cossolve/build/libcossolve.so'
try:
    libcossolve = ctypes.cdll.LoadLibrary(lib_path)
except:
    print ('Failed to load ' + lib_path)
    raise

# Properly set the return types for the C bindings
libcossolve.cossolve_StaticSolver_construct.restype = ctypes.c_voidp
libcossolve.cossolve_StaticSolver_delete.restype = None
libcossolve.cossolve_StaticSolver_addPointForce.restype = None
libcossolve.cossolve_StaticSolver_addDistributedForce.restype = None
libcossolve.cossolve_StaticSolver_clearForces.restype = None
libcossolve.cossolve_StaticSolver_addFixedConstraint.restype = None
libcossolve.cossolve_StaticSolver_getCoords.restype = ctypes.c_voidp
libcossolve.cossolve_StaticSolver_getStrains.restype = ctypes.c_voidp
libcossolve.cossolve_StaticSolver_getTwists.restype = ctypes.c_voidp
libcossolve.cossolve_StaticSolver_getFixedConstraintForces.restype = ctypes.c_voidp
libcossolve.cossolve_StaticSolver_getSystemMatrix.restype = ctypes.c_voidp
libcossolve.cossolve_StaticSolver_getSystemRows.restype = ctypes.c_int
libcossolve.cossolve_StaticSolver_getSystemCols.restype = ctypes.c_int
libcossolve.cossolve_StaticSolver_solve.restype = None
libcossolve.cossolve_StaticSolver_convergenceParameter.restype = ctypes.c_double

# Import objects so they can be accessed by anyone who includes this package
from .Solver import Solver
