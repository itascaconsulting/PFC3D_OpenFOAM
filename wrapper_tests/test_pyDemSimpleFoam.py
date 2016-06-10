import numpy as np
from pyDemFoam import pyDemSimpleFoam

solver = pyDemSimpleFoam()

assert  solver.nFaces() == 5616
assert solver.nNodes() == 2197
assert solver.nCells() == 1728
print solver.faces()
print solver.cell_faces()
print solver.face_centers()
print solver.phi()
print solver.phi().shape

x,y,z = solver.cell_centers().T
print 1
n=solver.n()
print 2
n[abs(y-0.005)<1e-9]=0.5
print (abs(y-0.005)<1e-9).sum()
solver.n(n)
print 3

solver.solve(0.0)
print 4

cc = solver.cell_centers()
flux = solver.cell_fluxes()
cv = solver.cell_volumes()

clist = range(120,130)

print solver.U()[:,1][clist]

for i in range(120,130):
    #print cc[i]
    print i
    print flux[i]
    print flux[i].sum()/cv[i]

assert type(solver.beta()) is np.ndarray
assert solver.beta().shape == (solver.nCells(),)
assert solver.ubar().shape == (solver.nCells(),3)
assert solver.U().shape == (solver.nCells(),3)
assert solver.n().shape == (solver.nCells(),)
assert solver.gradp().shape == (solver.nCells(),3)


print solver.flux_on_patch("topWall")
print solver.flux_on_patch("fixedWalls")
