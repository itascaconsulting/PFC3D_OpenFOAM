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
nmin, nmax = np.amin(nodes,axis=0), np.amax(nodes,axis=0)
diag = np.linalg.norm(nmin-nmax)
dmin, dmax = nmin -0.1*diag, nmax+0.1*diag

it.command("""
new
set deterministic on
set random 1
domain extent {} {} {} {} {} {}
set timestep max 1e-4
ball generate rad 0.00075 number 2160 box 0 0.04 0 0.06 0 0.0075
wall generate box 0 0.04 0 0.2 0 0.0075
ball ini dens 2000
cmat default model linear property kn 5000 ks 5000 fric 0.15 dp_nratio 0.1
set gravity 0 -9.81 0
def ball_height
  local max = 0
  loop foreach local b ball.list
    if ball.pos.y(b) > max then
      max = ball.pos.y(b)
    endif
  endloop
  ball_height = max
end
define pressure_drop
    pressure_drop = cfd_pressure_drop 
end
def fluid_time
  global fluid_time = mech.age
end
history add id 1 fish @fluid_time
history add id 2 fish @ball_height
plot clear
plot add ball
plot add wall transparency 70
cy 20000
ball ini damp 0.7
cy 20000
""".format(dmin[0], dmax[0],
           dmin[1], dmax[1],
           dmin[2], dmax[2]))
ca.create_mesh(nodes, elements)
it.command("""
config cfd
element cfd ini density {}
element cfd ini visc {}
cfd porosity poly
plot add cfdelement shape arrow colorby vectorattribute "velocity"
history add id 3 fish @pressure_drop
plot add hist 3 vs 1
""".format(fluid_density, fluid_viscosity))

element_volume = ca.volume()
dt = 0.0001

for i in range(100):
    it.command("solve time {}".format(dt))
    cfd_link.send_data(dt) # solve interval
    cfd_link.send_data(ca.porosity())
    cfd_link.send_data((ca.drag().T/element_volume).T/fluid_density)
    ca.set_pressure(cfd_link.read_data())
    ca.set_pressure_gradient(cfd_link.read_data())
    ca.set_velocity(cfd_link.read_data())
    it.fish.set("cfd_pressure_drop",(ca.pressure()[0] - ca.pressure()[1875]))

cfd_link.send_data(0.0) # solve interval
cfd_link.close()
del cfd_link

it.command("history write 1,2,3 file 'fluidized_bed_1.txt' truncate")
