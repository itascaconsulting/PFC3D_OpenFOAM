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

it.command("call ../fluidized_bed_4/make_ini.p3dat")
ca.create_mesh(nodes, elements)

it.command("""
config cfd
element cfd ini density {}
element cfd ini visc {}
cfd porosity poly
cfd buoy on
cfd interval 0
element cfd history id 3 zvelocity id 1
plot add hist 2 vs 1
""".format(fluid_density, fluid_viscosity))

element_volume = ca.volume()
oldu = ca.velocity()
old_force = (ca.drag().T/element_volume).T/fluid_density
old_porosity = ca.porosity()
r_factor = 0.5

smallest_size = 0.00095
time = 0.0
total_time = 4.0
max_dt = 0.0001

while time < total_time:
    dt = 100.0
    if (ca.velocity().max()>1.0e-7):
        dt = float(0.5*smallest_size/ca.velocity().max())
    if (max_dt < dt):
        dt = max_dt
    it.command("solve time {}".format(dt))
    it.command("cfd update")
    cfd_link.send_data(dt) # solve interval
    
    new_porosity = ca.porosity()
    cfd_link.send_data(old_porosity*(1-r_factor)+new_porosity*r_factor)
    old_porosity = new_porosity
    
    new_force = (ca.drag().T/element_volume).T/fluid_density
    cfd_link.send_data(old_force*(1-r_factor)+new_force*r_factor)
    old_force = new_force
    
    ca.set_pressure(cfd_link.read_data())
    ca.set_pressure_gradient(cfd_link.read_data())
    
    time += dt
    newu = cfd_link.read_data()
    ca.set_velocity(oldu*(1-r_factor)+newu*r_factor)
    oldu = newu

cfd_link.send_data(0.0) # solve interval
cfd_link.close()
del cfd_link

it.command("history write 1,2 file 'fluidized_bed_4.txt' truncate")
