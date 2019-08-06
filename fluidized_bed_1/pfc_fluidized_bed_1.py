import itasca as it
from itasca import cfdarray as ca
from itasca.util import p2pLinkServer

import numpy as np

cfd_link = p2pLinkServer()
cfd_link.start()

nodes = cfd_link.read_data()
elements = cfd_link.read_data()
fluid_density = cfd_link.read_data()
fluid_viscosity = cfd_link.read_data()
print fluid_density, fluid_viscosity
nmin, nmax = np.amin(nodes,axis=0), np.amax(nodes,axis=0)
diag = np.linalg.norm(nmin-nmax)
dmin, dmax = nmin -0.1*diag, nmax+0.1*diag
print dmin, dmax

it.command("""
new
domain extent {} {} {} {} {} {}
""".format(dmin[0], dmax[0],
           dmin[1], dmax[1],
           dmin[2], dmax[2]))
ca.create_mesh(nodes, elements)
it.command("""
config cfd
set timestep max 1e-5
element cfd ini density {}
element cfd ini visc {}
;cfd porosity poly
cfd buoy on
ball create rad 0.005 x 0.5 y 0.5 z 0.5
ball ini dens 2500
ball prop kn 1e2 ks 1e2 fric 0.25
set gravity 0 -9.81 0
def fluid_time
  global fluid_time = mech.age
end
history add id 1 fish @fluid_time
ball history id 2 yvelocity id 1
ball history id 3 yunbalforce id 1
ball cfd history id 4 yforce id 1
plot clear
plot add hist 2 vs 1
plot add cfdelement shape arrow colorby vectorattribute "velocity"
""".format(fluid_density, fluid_viscosity))

element_volume = ca.volume()
dt = 0.005

for i in range(100):
    it.command("solve age {}".format(it.mech_age()+dt))
    cfd_link.send_data(dt) # solve interval
    cfd_link.send_data(ca.porosity())
    cfd_link.send_data((ca.drag().T/element_volume).T/fluid_density)
    ca.set_pressure(cfd_link.read_data())
    ca.set_pressure_gradient(cfd_link.read_data())
    ca.set_velocity(cfd_link.read_data())

cfd_link.send_data(0.0) # solve interval
cfd_link.close()
del cfd_link

print "ball y velocity", it.ball.find(1).vel_y()
it.command("history write 1,2,3,4 file 'droptest1.txt' truncate")
