import numpy as np
from pyDemFoam import pyDemIcoFoam

solver = pyDemIcoFoam()

assert  solver.nFaces() == 5616
assert solver.nNodes() == 2197
assert solver.nCells() == 1728
print solver.faces()
print solver.cell_faces()
print solver.face_centers()
print solver.phi()
print solver.phi().shape


x,y,z = solver.cell_centers().T

n=solver.n()
n[abs(y-0.005)<1e-9]=0.5
print (abs(y-0.005)<1e-9).sum()
solver.n(n)

solver.solve(0.005)

cc = solver.cell_centers()
flux = solver.cell_fluxes()
cv = solver.cell_volumes()

clist = range(120,130)

print solver.U()[:,1][clist]

for i in range(120,130):
    #print cc[i]
    flux[i]
    flux[i].sum()*solver.dt()/cv[i]

print solver.flux_on_patch("topWall")
print solver.flux_on_patch("fixedWalls")
