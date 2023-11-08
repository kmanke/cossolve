import cossolve
import numpy

import matplotlib.pyplot as plt

solver = cossolve.Solver(4, 100e6, .3, 0, 2*.01**4/12, .01**4/12, .01**4/12, .01**2, 1, 0.005)
#F = -0.0001
#force = numpy.array([0, 0, F, 0, 0, 0])

#solver.add_force(1, force, False)

#x = numpy.linspace(0, 1, 100);
#y = F*x**2*(3-x)/(6*100e6*.01**4/12);

#fig = plt.figure()
#for it in range(0, 10):

#    solver.solve_strains()
    #solver.solve_coords()
#    g = solver.get_coords()

#    fig.clear()
#    plt.ylim([-0.5, 0.5])

   # plt.show(block = True)
#    fig.canvas.draw()
#    fig.canvas.flush_events()

#plt.plot(g[:,0,3], g[:,2,3], '-b')
#plt.plot(x, y, '--r')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('Comparison of Cosserat Rod vs. Bernoulli-Euler, F = {} N'.format(F))
#plt.show()
