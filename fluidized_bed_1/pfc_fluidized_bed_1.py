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

it.command("""
restore ../fluidized_bed_1/ini.sav
set timestep max 1e-4
set mechanical age 0.0
""")

ca.create_mesh(nodes, elements)

it.command("""
config cfd
element cfd ini density {}
element cfd ini visc {}
cfd porosity poly
history add id 3 fish @pressure_drop
plot clear
plot add ball
plot add axes
plot add wall transparency 70
plot add hist 3 vs 1
plot add cfdelement shape arrow colorby vectorattribute "velocity"
""".format(fluid_density, fluid_viscosity))

element_volume = ca.volume()
dt = 0.0001

for i in range(200):
    it.command("solve time {}".format(dt))
    cfd_link.send_data(dt) # solve interval
    cfd_link.send_data(ca.porosity())
    cfd_link.send_data((ca.drag().T/element_volume).T/fluid_density)
    ca.set_pressure(cfd_link.read_data())
    ca.set_pressure_gradient(cfd_link.read_data())
    ca.set_velocity(cfd_link.read_data())
    it.fish.set("cfd_pressure_drop",(ca.pressure()[12] - ca.pressure()[1887]))

cfd_link.send_data(0.0) # solve interval
cfd_link.close()
del cfd_link

#it.command("history write 1,2,3 file 'fluidized_bed_1.txt' truncate")
