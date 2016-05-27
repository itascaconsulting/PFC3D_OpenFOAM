import numpy as np
from pyDemSimpleFoam import pyDemSimpleFoam

solver = pyDemSimpleFoam()

# print solver.nFaces()
# print solver.nNodes()
# print solver.nCells()
# print solver.faces()
# print solver.cell_faces()
# print solver.face_centers()

# print solver.phi()
# print solver.phi().shape



x,y,z = solver.cell_centers().T

n=solver.n()
n[abs(y-0.005)<1e-9]=0.5
print (abs(y-0.005)<1e-9).sum()
solver.n(n)

solver.solve()

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

print solver.flux_on_patch("inletWalls")
print solver.flux_on_patch("outletWalls")
