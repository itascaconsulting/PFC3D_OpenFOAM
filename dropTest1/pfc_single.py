import itasca as it
it.command("python-reset-state false")
from itasca import cfdarray as ca
from pyDemFoam import pyDemIcoFoam
import numpy as np

solver = pyDemIcoFoam()
nodes = solver.nodes()
elements = solver.elements()
fluid_density = solver.rho()
fluid_viscosity = solver.mu()
print(fluid_density, fluid_viscosity)
nmin, nmax = np.amin(nodes,axis=0), np.amax(nodes,axis=0)
diag = np.linalg.norm(nmin-nmax)
dmin, dmax = nmin -0.1*diag, nmax+0.1*diag
print(dmin, dmax)

it.command("""
model new
MODEL LARGE-STRAIN on
domain extent {} {} {} {} {} {}
""".format(dmin[0], dmax[0],
           dmin[1], dmax[1],
           dmin[2], dmax[2]))
elements = elements.astype(np.longlong) # need to fix the c++ side, this is a work around for a platform type size issue.
ca.create_mesh(nodes, elements)

it.command("""
model config cfd
model mechanical timestep max 1e-5
element cfd ini density {}
element cfd ini visc {}
;cfd porosity poly
cfd buoy on
ball create rad 0.005 pos  0.5  0.5  0.5
ball ini dens 2500
ball prop "kn" 1e2 "ks" 1e2 "fric" 0.25
model gravity 0 0 -9.81
def fluid_time
  global fluid_time = mech.time.total
end
history fish id 1 fluid_time
ball history id 2 velocity-z id 1
ball history id 3 force-unbalanced-z id 1
;ball cfd history id 4 zforce id 1 ;; does not seem to work??
;plot clear
;plot add hist 2 vs 1
;plot add cfdelement shape arrow colorby vectorattribute "velocity"
""".format(fluid_density, fluid_viscosity))

element_volume = ca.volume()
dt = 0.005

for i in range(75):
    it.command("solve age {}".format(it.mech_age()+dt))
    cfd_link.send_data(dt) # solve interval
    solver.n(ca.porosity())
    solver.f((ca.drag().T/element_volume).T/fluid_density)
    solver.set_dt(dt)
    solver.solve(dt)

    ca.set_pressure(solver.p()*solver.rho())
    ca.set_pressure_gradient(solver.gradp()*solver.rho())
    ca.set_velocity(solver.U())

print("ball z velocity", it.ball.find(1).vel_z())
it.command("history export 1,2,3 file 'droptest1.txt' truncate")
it.command("model save 'final.sav'")
