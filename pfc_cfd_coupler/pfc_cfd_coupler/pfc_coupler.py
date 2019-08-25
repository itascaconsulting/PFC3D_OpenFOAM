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
        self.pressure = np.array([0]*self.nbElem)
        self.pressureMeasureCell1 = 0
        self.pressureMeasureCell2 = 1
        self.updatePressureDrop()
        self.max_dt = 0.005
        self.dt = self.max_dt
        self.time = 0.0
        self.cell_size = np.linalg.norm(self.elements_pos[0]-self.elements_pos[1])
        self.smallest_size = 1.0
        self.bandwidth = 2*self.cell_size
        
        nmin, nmax = np.amin(self.nodes,axis=0), np.amax(self.nodes,axis=0)
        diag = np.linalg.norm(nmin-nmax)
        dmin, dmax = nmin-0.1*diag, nmax+0.1*diag
        
        it.command("""
        new
        domain extent {} {} {} {} {} {}
        """.format(dmin[0], dmax[0], dmin[1], dmax[1], dmin[2], dmax[2]))

    def updatePressureDrop(self):
        it.fish.set("cfd_pressure_drop",(self.pressure[self.pressureMeasureCell1] - self.pressure[self.pressureMeasureCell2]))

    def updateWeights(self):
        bpos = ba.pos()
        btree = cKDTree(bpos,5)
        bmaps = btree.query_ball_tree(self.elements_tree,self.bandwidth)
        self.wmap=np.array([[0]*self.nbElem for x in bmaps],dtype='d')
        for ib in range(bpos.shape[0]):
            bp = bpos[ib]
            wlist = [0]*self.nbElem
            if len(bmaps[ib]):
                for ic in bmaps[ib]:
                    dv  = (vec((bp-self.elements_pos[ic]))).mag()
                    assert(dv<1)
                    wbc = self.kfunc(dv,self.bandwidth)
                    wlist[ic] = wbc
                self.wmap[ib] = wlist
                self.wmap[ib] /= self.wmap[ib].sum()
            else:
                d,iel = self.elements_tree.query(bp,k=1)
                self.wmap[ib,iel] = 1.0

    
    def kfunc(self,d,b):
        return math.exp(-(d/b)**2)
    
    def updatePorosity(self):
        # backward interpolations
        bvol = ba.mass_real() / ba.density()
        testdv = bvol - (self.wmap.T*bvol).sum(axis=0)
        assert(testdv.sum()<1e-20)
        
        # volume fraction and porosity in each cell
        evfracV = (self.wmap.T*bvol).sum(axis=1) 
        evfrac = evfracV / self.elements_vol
        self.elements_porosity = np.ones_like(evfrac) - evfrac

    def updateFluidDrag(self):
        bdrag = -1.0*self.balls_drag
        self.elements_drag = (np.einsum('ik,ij',self.wmap,bdrag)).T

    def updateBallsDrag(self):
        bvf   = np.einsum('ij,jk->ik',self.wmap,self.elements_vel)
        bvisc = np.einsum('ij,j...->i...'  ,self.wmap,self.elements_visc)
        bporo = np.einsum('ij,j...->i...'  ,self.wmap,self.elements_porosity)
        rel = bvf-ba.vel()
        self.balls_drag = np.full(rel.shape,0.0).T
        rho_f = self.fluid_density
        brad = ba.radius()
        brad2 = ba.radius()**2
        if rel.sum() != 0.0 :
            Reynolds = 2.0*rho_f*brad*np.linalg.norm(rel,axis=1)/ bvisc.T
            Cd = (0.63+4.8/np.sqrt(Reynolds))**2
            Chi = 3.7-0.65*np.exp(-(1.5-np.log10(Reynolds))**2/2.0)
            self.balls_drag = 0.5*rho_f*np.pi*brad2*Cd*rel.T*np.linalg.norm(rel,axis=1)*np.power(bporo,-1.0*Chi)
        self.balls_drag = self.balls_drag.T

    def updateForce(self):
        rho_f = self.fluid_density
        brad = np.array([ba.radius()]*3).T
        buoyancy = -4.0 / 3.0 * np.pi * brad**3 * rho_f * it.gravity()
        force = self.balls_drag + buoyancy
        ba.set_force_app(force)  

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

    def updateTimeStep(self):
        self.dt = 100.0
        if (self.elements_vel.max()>1.0e-7):
            self.dt = float(0.5*self.smallest_size/self.elements_vel.max())
        if (self.max_dt < self.dt):
            self.dt = self.max_dt

    def solve(self,total_time):
        self.time = 0.0
        self.balls_drag = np.full(ba.vel().shape,0.0)
        self.updateForce()
        while self.time < total_time:
            self.updateTimeStep()
            
            it.command("solve time {}".format(self.dt))
            
            self.updateWeights()
            self.updatePorosity()
            self.updateBallsDrag()
            self.updateFluidDrag()
            
            self.link.send_data(self.dt)
            self.link.send_data(self.elements_porosity)
            self.link.send_data((self.elements_drag.T/self.elements_vol).T/self.fluid_density)
            
            self.pressure = self.link.read_data()
            self.pressure_gradient = self.link.read_data()
            self.elements_vel = self.link.read_data()
            
            self.updatePressureDrop()
            self.updateBallsDrag()
            self.updateForce()
            
            self.time += self.dt

    def stopSolve(self):
        self.link.send_data(0.0)

    def close(self):
        self.link.close()
        del self.link    