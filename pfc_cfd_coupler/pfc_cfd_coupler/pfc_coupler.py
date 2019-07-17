import itasca as it
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
        self.elements_pos = self.link.read_data()
        self.elements_vol = self.link.read_data()
        self.fluid_density = self.link.read_data()
        self.fluid_viscosity = self.link.read_data()
        
        self.nbElem = self.elements.shape[0]
        self.elements_visc = np.array([self.fluid_viscosity]*self.nbElem)
        self.elements_tree = cKDTree(self.elements_pos)
        self.elements_vel = np.array([[0,0,0]]*self.nbElem)
        self.dt = 0.005
        
        nmin, nmax = np.amin(self.nodes,axis=0), np.amax(self.nodes,axis=0)
        diag = np.linalg.norm(nmin-nmax)
        dmin, dmax = nmin-0.1*diag, nmax+0.1*diag
        
        it.command("""
        new
        set timestep max 1e-5
        domain extent {} {} {} {} {} {}
        """.format(dmin[0], dmax[0], dmin[1], dmax[1], dmin[2], dmax[2]))

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
    
    def updatePorosityAndDrag(self):
        self.update_weights()
        
        # backward interpolations
        bvol = ba.mass_real() / ba.density()
        testdv = bvol - (self.wmap.T*bvol).sum(axis=0)
        assert(testdv.sum()<1e-20)
        
        # volume fraction and porosity in each cell
        evfracV = (self.wmap.T*bvol).sum(axis=1) 
        evfrac = evfracV / self.elements_vol
        self.elements_porosity = np.ones_like(evfrac) - evfrac
        
        bdrag = -1.0*ba.force_unbal()
        self.elements_drag = (np.einsum('ik,ij',self.wmap,bdrag)/self.elements_vol).T

    def updateForce(self):
        self.update_weights()
        
        #forward interpolations:
        
        self.testcv = (self.wmap*self.elements_vol).sum(axis=0)
        # testcv is equal to evol for each populated cell 
        
        #bvf is the fluid velocity for each ball
        bvf   = np.einsum('ij,jk->ik',self.wmap,self.elements_vel)
        bvisc = np.einsum('ij,j...->i...'  ,self.wmap,self.elements_visc)
        bporo = np.einsum('ij,j...->i...'  ,self.wmap,self.elements_porosity)
        
        self.bvf = bvf
        #ba.set_extra(1,bvf)
        rho_f = self.fluid_density
        brad = ba.radius()
        brad2 = ba.radius()**2
        
        buoyancy = np.zeros((brad2.shape[0],3))
        buoyancy[:,2] = -4.0 / 3.0 * np.pi * brad**3 * rho_f * it.gravity_z()
        
        ball_vel = ba.vel()
        rel = (bvf-ball_vel)
        force = np.full(ba.force_app().shape,0.0).T
        if rel.sum() != 0.0 :
            Reynolds = 2.0*rho_f*brad*np.linalg.norm(rel.T)/ bvisc.T
            Cd = (0.63+4.8/np.sqrt(Reynolds))**2
            Chi = 3.7-0.65*np.exp(-(1.5-np.log10(Reynolds))**2/2.0)
            force = 0.5*rho_f*np.pi*brad2*Cd*rel.T*np.linalg.norm(rel.T)*np.power(bporo,-1.0*Chi) + buoyancy.T
        
        ba.set_force_app(force.T)  

    #to see directions
    def plotFluidUnitVel(self):
        norms = np.sum(np.abs(self.elements_vel)**2,axis=-1)**(1./2)
        normsT = np.array([norms]).T
        arr = np.concatenate((self.elements_pos, self.elements_vel/normsT), 1)
        np.savetxt('vel.txt', arr, fmt='%1.6e', header="ITASCA VECTOR3D", comments='')
        it.command("vector import 'vel.txt'")

    #scaled by magnitude
    def plotFluidVel(self):
        arr = np.concatenate((self.elements_pos, self.elements_vel), 1)
        np.savetxt('vel.txt', arr, fmt='%1.6e', header="ITASCA VECTOR3D", comments='')
        it.command("vector import 'vel.txt'")

    def solve(self,nsteps):
        self.updatePorosityAndDrag()
        self.updateForce()
        for i in range(nsteps):
            it.command("solve age {}".format(it.mech_age()+self.dt))
            self.updatePorosityAndDrag()
            self.link.send_data(self.dt) # solve interval
            self.link.send_data(self.elements_porosity)
            self.link.send_data((self.elements_drag.T/self.elements_vol).T/self.fluid_density)
            self.pressure = self.link.read_data()
            self.pressure_gradient = self.link.read_data()
            self.elements_vel = self.link.read_data()
            self.updateForce()
        self.stopSolve()

    def stopSolve(self):
        self.link.send_data(0.0)

    def close(self):
        self.link.close()
        del self.link    