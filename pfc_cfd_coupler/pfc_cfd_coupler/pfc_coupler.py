import itasca as it
from itasca import cfdarray as ca
from itasca import ballarray as ba
from itasca.util import p2pLinkServer
import numpy as np
from scipy.spatial import cKDTree
from vec import vec
import math

class pfc_coupler(object):
    def __init__(self):
        self.link = p2pLinkServer()
        self.link.start()
        
        self.nodes = self.link.read_data()
        self.elements = self.link.read_data()
        self.nbElem = self.elements.shape[0]
        self.elements_pos = self.link.read_data()
        self.elements_vol = self.link.read_data()
        self.fluid_density = self.link.read_data()
        self.fluid_viscosity = self.link.read_data()
        self.elements_tree = cKDTree(self.elements_pos)
        
        #print fluid_density, fluid_viscosity
        nmin, nmax = np.amin(self.nodes,axis=0), np.amax(self.nodes,axis=0)
        diag = np.linalg.norm(nmin-nmax)
        dmin, dmax = nmin-0.1*diag, nmax+0.1*diag
        #print dmin, dmax
        
        it.command("""
        new
        domain extent {} {} {} {} {} {}
        """.format(dmin[0], dmax[0],
                   dmin[1], dmax[1],
                   dmin[2], dmax[2]))
        ca.create_mesh(self.nodes, self.elements)
        it.command("""
        config cfd
        set timestep max 1e-5
        element cfd ini density {}
        element cfd ini visc {}
        """.format(self.fluid_density, self.fluid_viscosity))

    def update_weights(self):
        bandwidth = 0.15
        bpos = ba.pos()
        btree = cKDTree(bpos,5)
        bmaps = btree.query_ball_tree(self.elements_tree,bandwidth)
        self.wmap=np.array([[None]*self.nbElem for x in bmaps],dtype='d')
        for ib in range(bpos.shape[0]):
          bp = bpos[ib]
          wlist = [0]*self.nbElem
          for ic in bmaps[ib]:
            dv  = (vec((bp-self.elements_pos[ic]))).mag()
            assert(dv<1)
            wbc = self.kfunc(dv,bandwidth)
            wlist[ic] = wbc
          self.wmap[ib] = wlist
        self.wmap /= self.wmap.sum(axis=1, keepdims=True)
    
    def kfunc(self,d,b):
        return math.exp(-(d/b)**2)
    
    def updatePorosity(self):
        #update porosity
        print "Updating Porosity..."
        
        self.update_weights()
        
        # backward interpolations
        bvol = ba.mass_real() / ba.density()
        testdv = bvol - (self.wmap.T*bvol).sum(axis=0)
        assert(testdv.sum()<1e-20)
        
        # volume fraction and porosity in each cell
        evfracV = (self.wmap.T*bvol).sum(axis=1) 
        evfrac = evfracV / self.elements_vol
        self.elements_porosity = np.ones_like(evfrac) - evfrac

    def solve(self):
        element_volume = ca.volume()
        dt = 0.005
        
        for i in range(100):
            it.command("solve age {}".format(it.mech_age()+dt))
            self.updatePorosity()
            print "sending solve time"
            self.link.send_data(dt) # solve interval
            self.link.send_data(self.elements_porosity)
            self.link.send_data((ca.drag().T/element_volume).T/self.fluid_density)
            print " cfd solve started"
            ca.set_pressure(self.link.read_data())
            ca.set_pressure_gradient(self.link.read_data())
            ca.set_velocity(self.link.read_data())
            print " cfd solve ended"
        self.link.send_data(0.0) # solve interval