import cossolve
import numpy
import sys

import matplotlib.pyplot as plt

numpy.set_printoptions(threshold=sys.maxsize, precision=10, linewidth=300, suppress=False)

num_nodes = 50
num_its = 5
solver = cossolve.Solver(num_nodes, 100e6, .3, 0, 2*.01**4/12, .01**4/12, .01**4/12, .01**2, 1, 0.005)

# Add a fixed constraint at the left end
g = numpy.eye(4, order='F')
g[0, 3] = 1-1.0/49
g[2, 3] = -.1
#solver.add_fixed_constraint(49, g)

print(2*numpy.pi*100e6*.01**4/12/1)
F = -.5236
force = numpy.array([0, 0, 0, 0, F, 0])
solver.add_point_force(1, force, True)

for i in range(0, num_its):
    solver.solve()

strains = solver.get_strains()
twists = solver.get_twists()
#reactions = solver.get_fixed_constraint_forces()
sysmat = solver.get_system_matrix()
coords = solver.get_coords()

x = numpy.arange(0, 1, 1/49);
y = F*x**2*(3-x)/(6*100e6*.01**4/12);

fig = plt.figure()

plt.plot(coords[:,0,3], coords[:,2,3], '-b')
#plt.plot(x, y, '--r')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Model Validation - Euler-Bernoulli'.format(F))
plt.show()
