# Python interface for Solver class.
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
import numpy

from cossolve import libcossolve

class Solver:
    # Constructor: create the C++ object and store a handle
    def __init__(self, num_nodes, tensile_modulus, poisson_ratio, shear_modulus,
                 moment_x, moment_y, moment_z, area, length, linear_density):
        self.handle = ctypes.c_voidp(libcossolve.cossolve_createSolver(
            ctypes.c_int(num_nodes), ctypes.c_double(tensile_modulus),
            ctypes.c_double(poisson_ratio), ctypes.c_double(shear_modulus),
            ctypes.c_double(moment_x), ctypes.c_double(moment_y),
            ctypes.c_double(moment_z), ctypes.c_double(area), ctypes.c_double(length),
            ctypes.c_double(linear_density)
            ))

    # Destructor: deallocate the solver
    def __del__(self):
        libcossolve.cossolve_deleteSolver(self.handle)

    # Returns the number of nodes in this solver
    @property
    def num_nodes(self):
        return ctypes.c_int(libcossolve.cossolve_getNodeCount(self.handle)).value
        
    # Returns the strain vector as a numpy array.
    def get_strains(self):
        strains = numpy.empty((self.num_nodes * 6, 1))
        libcossolve.cossolve_getStrains(self.handle, ctypes.c_voidp(strains.ctypes.data))
        return strains

    # Returns the coordinate transformations as a numpy array.
    def get_coords(self):
        coords = numpy.empty((self.num_nodes, 4, 4))
        libcossolve.cossolve_getCoords(self.handle, ctypes.c_voidp(coords.ctypes.data))

        # Set the proper data order
        coords.strides = (128, 8, 32)
        return coords

    # Adds the specified force to the solver
    def add_force(self, s, force, body_frame):
        libcossolve.cossolve_Solver_addForce(self.handle, ctypes.c_double(s),
                                             ctypes.c_voidp(force.ctypes.data),
                                             ctypes.c_bool(body_frame))
        return

    def solve_strains(self):
        libcossolve.cossolve_Solver_solveStrains(self.handle)
        return

    def solve_coords(self):
        libcossolve.cossolve_Solver_solveCoords(self.handle)
        return
