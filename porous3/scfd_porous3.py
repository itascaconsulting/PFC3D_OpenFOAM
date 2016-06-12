import numpy as np
from pyDemFoam import pyDemSimpleFoam2, pyDemSimpleFoam

solver = pyDemSimpleFoam()

nc = solver.nCells()
x,y,z = solver.cell_centers().T

solver.n(np.ones(nc))
solver.beta(np.zeros(nc))
solver.ubar(np.zeros((nc,3)))

solver.solve(0.0)
solver.solve(0.0)
solver.solve(0.0)
solver.solve(0.0)

print solver.U()
