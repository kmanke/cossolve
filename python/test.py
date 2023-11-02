import cossolve
import numpy

import matplotlib.pyplot as plt

solver = cossolve.Solver(50, 100e6, .3, 0, 2*.01**4/12, .01**4/12, .01**4/12, .01**2, 1, 0.1)
#force = numpy.array([0, 0, -0.2, 0, 0, 0])

#solver.add_force(1, force, True)

while(True):

    solver.solve_strains()
    #solver.solve_coords()
    g = solver.get_coords()
    plt.plot(g[:,0,3], g[:,2,3])
    plt.ylim([-0.5, 0.5])

    plt.show()
