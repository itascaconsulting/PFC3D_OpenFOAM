import itasca as it
from itasca import cfdarray as ca
from itasca.util import p2pLinkServer

import numpy as np

with p2pLinkServer() as cfd_link:
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
    element cfd ini density {}
    element cfd ini visc {}
    cfd porosity poly
    cfd buoy on
    cfd relax 1
    cfd update
    call particles.p3dat
    """.format(fluid_density, fluid_viscosity))


    element_volume = ca.volume()
    dt = 0.005
    oldu = ca.velocity()
    r_factor = 1.0
    for i in range(20):
        it.command("solve age {}".format(it.mech_age()+dt))
        print "sending solve time"
        cfd_link.send_data(dt) # solve interval
        cfd_link.send_data(ca.porosity())
        new_force = (ca.drag().T/element_volume).T/fluid_density
        cfd_link.send_data(new_force)
        print " cfd solve started"
        ca.set_pressure(cfd_link.read_data())
        ca.set_pressure_gradient(cfd_link.read_data())
        newu = cfd_link.read_data()
        ca.set_velocity(oldu*(1-r_factor)+newu*r_factor)
        oldu = newu
        print " cfd solve ended"

    cfd_link.send_data(0.0)

    print "ball z velocity", it.ball.find(1).vel_z()
it.command('plot bitmap filename "porous1.png"')
