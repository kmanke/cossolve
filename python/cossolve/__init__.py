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
libcossolve.cossolve_createSolver.restype = ctypes.c_voidp
libcossolve.cossolve_deleteSolver.restype = None
libcossolve.cossolve_getStrains.resType = None
libcossolve.cossolve_getCoords.resType = None
libcossolve.cossolve_getNodeCount.resType = ctypes.c_int
libcossolve.cossolve_Solver_addForce.resType = None
libcossolve.cossolve_Solver_solveStrains.resType = None
libcossolve.cossolve_Solver_solveCoords.resType = None

# Import objects so they can be accessed by anyone who includes this package
from .Solver import Solver
