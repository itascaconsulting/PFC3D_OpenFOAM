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
    set timestep max 1e-5
    element cfd ini density {}
    element cfd ini visc {}
    cfd porosity poly
    cfd buoy on
    ball create rad 0.005 x 0.5 y 0.5 z 0.5
    ball ini dens 2500
    ball prop kn 1e2 ks 1e2 fric 0.25
    set gravity 0 0 -9.81
    def fluid_time
      global fluid_time = mech.age
    end
    ball history id 1 zvelocity id 1
    history add id 2 fish @fluid_time
    plot clear
    plot add hist 1 vs 2
    plot add cfdelement shape arrow colorby vectorattribute "velocity"
    """.format(fluid_density, fluid_viscosity))

    element_volume = ca.volume()
    dt = 0.005

    for i in range(100):
        it.command("solve age {}".format(it.mech_age()+dt))
        print "sending solve time"
        cfd_link.send_data(dt) # solve interval
        cfd_link.send_data(ca.porosity())
        cfd_link.send_data((ca.drag().T/element_volume).T/fluid_density)
        print " cfd solve started"
        ca.set_pressure(cfd_link.read_data())
        ca.set_pressure_gradient(cfd_link.read_data())
        ca.set_velocity(cfd_link.read_data())
        print " cfd solve ended"

    cfd_link.send_data(0.0) # solve interval


    print "ball z velocity", it.ball.find(1).vel_z()
