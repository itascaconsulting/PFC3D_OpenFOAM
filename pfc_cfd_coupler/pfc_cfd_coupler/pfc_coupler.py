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
        
        nmin, nmax = np.amin(self.nodes,axis=0), np.amax(self.nodes,axis=0)
        diag = np.linalg.norm(nmin-nmax)
        dmin, dmax = nmin-0.1*diag, nmax+0.1*diag
        self.corner1 = nmin
        self.corner2 = nmax
        self.grid_sizes = np.array([20,10,10])
        self.cell_sizes = (self.corner2-self.corner1)/self.grid_sizes
        self.cell_size = np.linalg.norm(self.cell_sizes)
        self.smallest_size = 1.0
        self.bandwidth = 2*self.cell_size
        
        it.command("""
        new
        domain extent {} {} {} {} {} {}
        """.format(dmin[0], dmax[0], dmin[1], dmax[1], dmin[2], dmax[2]))

    def updatePressureDrop(self):
        it.fish.set("cfd_pressure_drop",(self.pressure[self.pressureMeasureCell1] - self.pressure[self.pressureMeasureCell2]))

    def updateWeights(self):
        bpos = ba.pos()
        nb = bpos.shape[0]
        nmin = np.array([self.corner1]*nb)
        nmax = np.array([self.corner2]*nb)
        dmin = np.concatenate((bpos-nmin,nmax-bpos),1)
        dminx = dmin[:,0:6:3]
        dminy = dmin[:,1:6:3]
        dminz = dmin[:,2:6:3]
        self.d = np.vstack((dminx.min(axis=1),dminy.min(axis=1),dminz.min(axis=1))).T
        btree = cKDTree(bpos,5)
        bmaps = btree.query_ball_tree(self.elements_tree,self.bandwidth)
        self.wmap=np.array([[0]*self.nbElem for x in bmaps],dtype='d')
        for ib in range(nb):
            bp = bpos[ib]
            wlist = [0]*self.nbElem
            if len(bmaps[ib]):
                for ic in bmaps[ib]:
                    dv  = (vec((bp-self.elements_pos[ic]))).mag()
                    assert(dv<1)
                    wbc = self.kfunc(dv,self.bandwidth,self.d[ib])
                    wlist[ic] = wbc
                self.wmap[ib] = wlist
                self.wmap[ib] /= self.wmap[ib].sum()
            else:
                d,iel = self.elements_tree.query(bp,k=1)
                self.wmap[ib,iel] = 1.0
    
    def kfunc(self,d,b,a):
        x = d/b
        s = self.smallest_size
        pow1 = self.pow1
        pow2 = self.pow2
        pow3 = self.pow3
        if x<1:
            if   a[0]>=s and a[1]>=s and a[2]>=s:
                return (1-x**2)**4
            elif a[0]<s  and a[1]>=s and a[2]>=s:
                return (1-x**2)**(-4./s*a[0]+pow1)
            elif a[0]>=s and a[1]<s  and a[2]>=s:
                return (1-x**2)**(-4./s*a[1]+pow1)
            elif a[0]>=s and a[1]>=s and a[2]<s:
                return (1-x**2)**(-4./s*a[2]+pow1)
            elif a[0]<s  and a[1]<s  and a[2]>=s:
                p = -4./s*a[0]+pow2
                return (1-x**2)**(-p/s*a[1]+pow2)
            elif a[0]<s  and a[1]>=s and a[2]<s:
                p = -4./s*a[2]+pow2
                return (1-x**2)**(-p/s*a[0]+pow2)
            elif a[0]>=s and a[1]<s  and a[2]<s:
                p = -4./s*a[1]+pow2
                return (1-x**2)**(-p/s*a[2]+pow2)
            else:
                p = -4./s*a[0]+pow3
                q = -p /s*a[1]+pow3
                return (1-x**2)**(-q/s*a[2]+pow3)
        else:
            return 0
    
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
        bporo = np.einsum('ij,j...->i...'  ,self.wmap,self.elements_porosity)
        rel = bvf-ba.vel()
        self.balls_drag = np.full(rel.shape,0.0).T
        rho_f = self.fluid_density
        brad = ba.radius()
        brad2 = ba.radius()**2
        if rel.sum() != 0.0 :
            Reynolds = 2.0*rho_f*brad*np.linalg.norm(rel,axis=1)/ self.fluid_viscosity
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

    def plotPorosity(self):
        f = open("elements.geom","w")
        f.write("ITASCA GEOMETRY3D\nNODES\n")
        for i in range(1,self.nbElem+1):
            f.write(str(i))
            f.write(' ')
            for j in self.elements_pos[i-1]:
                f.write(str(j))
                f.write(' ')
            f.write('EXTRA 1 ')
            f.write(str(self.elements_porosity[i-1]))
            f.write(' EXTRA 2 ')
            f.write(str(np.linalg.norm(self.elements_drag[i-1])))
            f.write('\n')
        f.close()
        it.command("geometry import 'elements.geom'")

    def updateTimeStep(self):
        self.dt = 100.0
        if (self.elements_vel.max()>1.0e-7):
            self.dt = float(0.5*self.smallest_size/self.elements_vel.max())
        if (self.max_dt < self.dt):
            self.dt = self.max_dt

    def initialize(self):
        self.time = 0.0
        self.balls_drag = np.full(ba.vel().shape,0.0)
        self.updateForce()

    def solve(self,solve_time):
        time = self.time
        while self.time < time + solve_time:
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