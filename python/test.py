import cossolve
import numpy
import sys

import matplotlib.pyplot as plt

numpy.set_printoptions(threshold=sys.maxsize, precision=1, linewidth=300, suppress=False)

num_nodes = 100
num_its = 100

solver = cossolve.Solver(num_nodes, 100e6, .3, 0, 2*.01**4/12, .01**4/12, .01**4/12, .01**2, 1, 0.005)

# Add a fixed constraint at the left end
#g = numpy.eye(4, order='F')
#g[0, 3] = .99
#g[2, 3] = 0
#solver.add_fixed_constraint(99, g)

F = -1
force = numpy.array([0, 0, F*.5, 0, 0, 0])

solver.add_point_force(1, force, True)

for i in range(0, num_its):
    solver.solve()

strains = solver.get_strains()
twists = solver.get_twists()
#reactions = solver.get_fixed_constraint_forces()
sysmat = solver.get_system_matrix()
coords = solver.get_coords()

#print(reactions)

x = numpy.linspace(0, 1, 100);
y = F*x**2*(3-x)/(6*100e6*.01**4/12);

fig = plt.figure()

plt.plot(coords[:,0,3], coords[:,2,3], '-b')
#plt.plot(x, y, '--r')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Comparison of Cosserat Rod vs. Bernoulli-Euler, F = {} N'.format(F))
plt.show()
