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
        self.handle = ctypes.c_voidp(libcossolve.cossolve_StaticSolver_construct(
            ctypes.c_int(num_nodes), ctypes.c_double(tensile_modulus),
            ctypes.c_double(poisson_ratio), ctypes.c_double(shear_modulus),
            ctypes.c_double(moment_x), ctypes.c_double(moment_y),
            ctypes.c_double(moment_z), ctypes.c_double(area), ctypes.c_double(length),
            ctypes.c_double(linear_density)
            ))
        self.num_nodes = num_nodes

    # Destructor: deallocate the solver
    def __del__(self):
        libcossolve.cossolve_StaticSolver_delete(self.handle)

    # Returns the coordinate transformations as a numpy array.
    def get_coords(self):
        coords = numpy.empty((self.num_nodes, 4, 4))
        libcossolve.cossolve_StaticSolver_getCoords(self.handle, ctypes.c_voidp(coords.ctypes.data))

        # Set the proper data order
        coords.strides = (128, 8, 32)
        return coords
        
    # Returns the strain vector as a numpy array.
    def get_strains(self):
        strains = numpy.empty((self.num_nodes * 6))
        libcossolve.cossolve_StaticSolver_getStrains(self.handle, ctypes.c_voidp(strains.ctypes.data))
        return strains

    # Returns the twist vector as a numpy array.
    def get_twists(self):
        twists = numpy.empty((self.num_nodes * 6))
        libcossolve.cossolve_StaticSolver_getTwists(self.handle, ctypes.c_voidp(twists.ctypes.data))
        return twists

    def get_fixed_constraint_forces(self):
        forces = numpy.empty((6))
        libcossolve.cossolve_StaticSolver_getFixedConstraintForces(self.handle, ctypes.c_voidp(forces.ctypes.data))
        return forces
        
    def get_system_matrix(self):
        mat = numpy.empty((self.get_system_rows(), self.get_system_cols()), order='F')
        libcossolve.cossolve_StaticSolver_getSystemMatrix(self.handle, ctypes.c_voidp(mat.ctypes.data))
        return mat

    def get_system_rows(self):
        return libcossolve.cossolve_StaticSolver_getSystemRows(self.handle)

    def get_system_cols(self):
        return libcossolve.cossolve_StaticSolver_getSystemCols(self.handle)
    
    # Adds the specified force to the solver
    def add_point_force(self, s, force, body_frame):
        libcossolve.cossolve_StaticSolver_addPointForce(self.handle, ctypes.c_double(s),
                                                        ctypes.c_voidp(force.ctypes.data),
                                                        ctypes.c_bool(body_frame))
        return

    def add_distributed_force(self, s1, s2, force, body_frame):
        libcossolve.cossolve_StaticSolver_addDistributedForce(self.handle, ctypes.c_double(s1),
                                                              ctypes.c_double(s2),
                                                              ctypes.c_voidp(force.ctypes.data),
                                                              ctypes.c_bool(body_frame))
        return

    def clear_forces(self):
        libcossolve.cossolve_StaticSolver_clearForces(self.handle)
        return
    
    def add_fixed_constraint(self, node, g):
        libcossolve.cossolve_StaticSolver_addFixedConstraint(self.handle, ctypes.c_int(node),
                                                             ctypes.c_voidp(g.ctypes.data))
        return

    def solve(self):
        libcossolve.cossolve_StaticSolver_solve(self.handle)
        return
