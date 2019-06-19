import itasca as it
from itasca import cfdarray as ca
from itasca.util import p2pLinkServer
import numpy as np

class pfc_coupler(object):
    def __init__(self):
        self.link = p2pLinkServer()
        self.link.start()
        
        nodes = self.link.read_data()
        elements = self.link.read_data()
        self.fluid_density = self.link.read_data()
        fluid_viscosity = self.link.read_data()
        #print fluid_density, fluid_viscosity
        nmin, nmax = np.amin(nodes,axis=0), np.amax(nodes,axis=0)
        diag = np.linalg.norm(nmin-nmax)
        dmin, dmax = nmin-0.1*diag, nmax+0.1*diag
        #print dmin, dmax
        
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
        """.format(self.fluid_density, fluid_viscosity))

    def solve(self):
        element_volume = ca.volume()
        dt = 0.005
        
        for i in range(100):
            it.command("solve age {}".format(it.mech_age()+dt))
            print "sending solve time"
            self.link.send_data(dt) # solve interval
            self.link.send_data(ca.porosity())
            self.link.send_data((ca.drag().T/element_volume).T/self.fluid_density)
            print " cfd solve started"
            ca.set_pressure(self.link.read_data())
            ca.set_pressure_gradient(self.link.read_data())
            ca.set_velocity(self.link.read_data())
            print " cfd solve ended"
        self.link.send_data(0.0) # solve interval