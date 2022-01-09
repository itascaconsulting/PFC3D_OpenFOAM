import itasca
import math
import numpy
from scipy.spatial import cKDTree
from itasca import ballarray as ba
from vec import vec
import customsocket
import os

class ScTetraCoupler(object):
  def __init__(self,inistate):
    """restore the initial PFC state""" 
    itasca.command("restore {}".format(inistate))
    itasca.command("set mechanical age 0.0")
    itasca.command("ball extra 1 [vector(0.0,0.0,0.0)]")
    itasca.command("""
    define set_geometry_extra(gs_id,ex_idx,ex_val)
      local gs = geom.set.find(gs_id)
      loop foreach local p geom.poly.list(gs)
        geom.poly.extra(p,ex_idx) = ex_val
      endloop
    end
    """)
    pass

  def execute(self,tcp_id='127.0.0.1', tcp_port=3333):
    self.initialize(tcp_id=tcp_id,tcp_port=tcp_port)
    self.couple()
    self.terminate()

  def update_weights(self):
    bandwidth = 10.0e-3
    bpos = ba.pos()
    btree = cKDTree(bpos,5)
    bmaps = btree.query_ball_tree(self.elements_tree,bandwidth)
    self.wmap=numpy.array([[None]*self.nbElem for x in bmaps],dtype='d')
    for ib,clist in numpy.ndenumerate(bmaps):
      bp = bpos[ib]
      wlist = [0]*self.nbElem
      for ic in clist:
        dv  = (vec((bp-self.elements_pos[ic]))).mag()
        assert(dv<1)
        wbc = self.kfunc(dv,bandwidth)
        wlist[ic] = wbc
      self.wmap[ib] = wlist
    self.wmap /= self.wmap.sum(axis=1, keepdims=True)
    
  def kfunc(self,d,b):
    return math.exp(-(d/b)**2)
  
  def updatePorosityAndDrag(self):
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
    self.elements_porosity = numpy.ones_like(evfrac) - evfrac
    
    # drag force on volume element
    bdrag = -1.0*ba.force_app()
    self.elements_drag = (numpy.einsum('ik,ij',self.wmap,bdrag)/self.elements_vol).T
    
  def sendPorosityAndDrag(self):  
    print "Sending porosity to SCTetra..."
    for ie in range(0,self.nbElem):
      self.client.send_data(float(self.elements_porosity[ie]))
    print "Porosity sent to SCTetra !"
    print "Sending drag to SCTetra..."
    for ie in range(0,self.nbElem):
      for i in range(0,itasca.dim()):
        self.client.send_data(-0.001*float(self.elements_drag[ie][i]))
    print "Drag sent to SCTetra !"
      
  def receiveCFDTime(self):
    self.cfdDT   = self.client.read_data()
    self.cfdTime = self.client.read_data()

  def receiveCFDData(self):
    # UPDATE FIELD
    # receive number of elements
    print "Waiting for updated CFD data..."
    nbElem = self.client.read_data()
    print "Number of elements: ", nbElem
    # SE: what if nbElem is not consistent with self.nbElem ? 
    self.elements_field = {}
    for n in range(0,nbElem):
      u = self.client.read_data()
      v = self.client.read_data()
      w = self.client.read_data()
      p = self.client.read_data()
      dpx = self.client.read_data()
      dpy = self.client.read_data()
      dpz = self.client.read_data()
      d = self.client.read_data()
      vi = self.client.read_data()
      self.elements_field[n] = [u,v,w,p,dpx,dpy,dpz,d,vi]
        
    # Update data
    #store elements velocity into a numpy array
    a = numpy.array(self.elements_field.values())
    # remove last column
    self.elements_vel = numpy.delete(a,(3,4,5,6,7,8),1)
    self.elements_visc = numpy.delete(a,(0,1,2,3,4,5,6,7),1)
    #update wmap
    self.update_weights()
    
    # update pressure drop:
    itasca.fish.set("cfd_pressure_drop",(self.elements_field[self.gpid1][3] - self.elements_field[self.gpid2][3]))
    
    #forward interpolations:

    self.testcv = (self.wmap*self.elements_vol).sum(axis=0)
    # testcv is equal to evol for each populated cell 
    
    #bvf is the fluid velocity for each ball
    bvf   = numpy.einsum('ij,jk->ik',self.wmap,self.elements_vel)
    bvisc = numpy.einsum('ij,j...->i...'  ,self.wmap,self.elements_visc)
    bporo = numpy.einsum('ij,j...->i...'  ,self.wmap,self.elements_porosity)

    self.bvf = bvf
    ba.set_extra(1,bvf)
    #Cd=0.4
    rho_f = 1200.0
    brad = ba.radius()
    brad2 = ba.radius()**2

    #buoyancy = numpy.zeros((brad2.shape[0],3))
    #buoyancy[:,1] = -4.0 / 3.0 * numpy.pi * ba.radius()**3 * rho_f * itasca.gravity_y()
    
    ball_vel = ba.vel()
    rel = (bvf-ball_vel)
    force = numpy.full(ba.force_app().shape,0.0).T
    if rel.sum() != 0.0 :
      Reynolds = 2.0*rho_f*brad*abs(rel.T)/ bvisc.T
      Cd = (0.63+4.8/numpy.sqrt(Reynolds))**2
      Chi = 3.7-0.65*numpy.exp(-(1.5-numpy.log10(Reynolds))**2/2.0)
      force = 0.5*rho_f*numpy.pi*brad2*Cd*rel.T*abs(rel.T)*numpy.power(bporo,-1.0*Chi)  #+ buoyancy.T
    
    ba.set_force_app(force.T)    
    
    print "Received CFD data..."
    #send communication status (success)
    self.client.send_data(1)    

  def update_geom_PressureAndVel(self):
    itasca.command("set echo off")
    for ie in range(0,self.nbElem):
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,5,float(self.elements_field[ie][0])))  # x-vel
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,6,float(self.elements_field[ie][1])))  # y-vel
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,7,float(self.elements_field[ie][2])))  # z-vel
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,8,float(self.elements_field[ie][3])))  # pressure
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,9,float(self.elements_field[ie][4])))   # dpx
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,10,float(self.elements_field[ie][5])))  # dpy
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,11,float(self.elements_field[ie][6])))  # dpz
    itasca.command("set echo on")

  def update_geom_PoroAndDrag(self):
    #update geometry set extra
    itasca.command("set echo off")
    for ie in range(0,self.nbElem):
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,1,float(self.elements_porosity[ie])))
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,2,float(self.elements_drag[ie][0])))
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,3,float(self.elements_drag[ie][1])))
      itasca.command("@set_geometry_extra({},{},{})".format(ie+1,4,float(self.elements_drag[ie][2])))
    itasca.command("set echo on")

  def couple(self):
    # update porosity and drag  
    self.updatePorosityAndDrag()
    # receive initial CFD data
    self.receiveCFDData()
    # receive initial CFD time
    self.receiveCFDTime()
    self.NIter = 0
    while self.cfdTime >= 0.0:
      # solve to catch up CFD time
      # if CFD Time remains 0 ScTetra is 
      # in stationnary mode... How do we handle this ?
      if self.cfdTime > 0.0:
        # attempt to fix the timestep to a reasonable value
        #itasca.command("set timestep fix {}".format(self.cfdDT/100.0))
        self.NIter += 1
        print " --- CFD cycle {} --- ".format(self.NIter)
        itasca.command("solve time {} exact".format(self.cfdDT))
      # update porosity and drag  
      self.updatePorosityAndDrag()
      # send porosiy and drag
      self.sendPorosityAndDrag()
      # receive updated CFD data
      self.receiveCFDData()
      # receive updated CFD time
      self.receiveCFDTime()
    
    # terminate if cfdtime is negative
    print "received negative CFD time of {}".format(self.cfdTime)
    self.terminate()
    
  def terminate(self):
    self.client.close()
  
  def initialize(self,tcp_id='127.0.0.1', tcp_port=3333):
    # open channel (client mode)
    self.client = customsocket.CustomSocketClient(tcp_id=tcp_id,tcp_port=tcp_port)
    self.client.start()
    # receive mesh from SCTetra
    # receive number of nodes
    self.nbNodes = self.client.read_data()
    print "Number of nodes: ",self.nbNodes
    
    # receive node information
    self.nodes = {}
    for n in range(0,self.nbNodes):
      i = self.client.read_data() # node id
      x = self.client.read_data() # node x 
      y = self.client.read_data() # node y
      z = self.client.read_data() # node z
      self.nodes[i] = [x,y,z]
      #print "coord node {} : {} {} {}".format(n,x,y,z)

    # receive connectivity information
    self.nbElem = self.client.read_data()
    print "Number of elements: ",self.nbElem
    self.elements_nodes = {}
    for n in range(0,self.nbElem):
      i  = self.client.read_data()  # element ID (should be n)
      nn = self.client.read_data()  # total number of nodes
      idx = []
      for j in range(0,nn):
        nj = self.client.read_data()  # index of node j
        idx.append(nj)
      self.elements_nodes[i] = idx
    #  print "connection {} : {}".format(n,ic)
    
    # SE: mesh creation takes a while !
    self.create_mesh()
    
    #SE: find index of polygon closest to a specified position in all sets
    itasca.command("""
    define poly_near(v)
      local ret = null
      local d = 1e20
      loop foreach local gs geom.set.list
        local gp = geom.poly.near(gs,v)
        if gp # null then
          local dist = math.mag(geom.poly.pos(gp)-v)
          if dist < d then
            ret = gp
            d = dist
          endif
        endif
      endloop
      poly_near = geom.poly.id(ret)
    end
    [gp1 = poly_near(vector(2.0e-2,2.0e-3,3.75e-3))]
    [gp2 = poly_near(vector(2.0e-2,15.0e-2,3.75e-3))]
    
    define pressure_drop
      pressure_drop = cfd_pressure_drop 
    end
        
    history id 1 mechanical age
    history id 2 @pressure_drop
    """)
    
    self.gpid1 = itasca.fish.get("gp1") - 1
    self.gpid2 = itasca.fish.get("gp2") - 1
    
    # UPDATE MESHELEM
    # receive number of elements
    self.nbElem = self.client.read_data()
    print "Number of elements: ",self.nbElem
    self.elements_geom = {}
    for n in range(0,self.nbElem):
      x = self.client.read_data()
      y = self.client.read_data()
      z = self.client.read_data()
      v = self.client.read_data()
      self.elements_geom[n] = [x,y,z,v]
      #print "element {} : x = {} y = {} z = {} - vol = {}".format(n,x,y,z,v)
    
    #store elements centroids into a numpy array
    a = numpy.array(self.elements_geom.values())
    # remove last column
    self.elements_pos = numpy.delete(a,-1,1)
    # remove first 3 columns
    self.elements_vol = numpy.delete(a,(0,1,2),1).reshape(-1)
    
    # create a cKDTree for lookup
    self.elements_tree = cKDTree(self.elements_pos)
    
  def create_mesh(self):
    # write command structure into a temporary file
    file = open("~tmp.txt", "w")
    file.write("define create_mesh\n")
    nlines = 1
    for ie in range(0,len(self.elements_nodes)):
      s = "poly_{}".format(ie+1)
      file.write("  gset = geom.set.find({})\n".format(ie+1))
      file.write("  if gset # null then\n")
      file.write("    geom.set.delete(gset)\n")
      file.write("  endif\n")
      file.write("  gset = geom.set.create(\"{}\",{})\n".format(s,ie+1))
      nlines += 5
      l = self.elements_nodes[ie]
      if len(l) == 4:
        # tetra
        x0 = self.nodes[l[0]][0]
        y0 = self.nodes[l[0]][1]
        z0 = self.nodes[l[0]][2]
        x1 = self.nodes[l[1]][0]
        y1 = self.nodes[l[1]][1]
        z1 = self.nodes[l[1]][2]
        x2 = self.nodes[l[2]][0]
        y2 = self.nodes[l[2]][1]
        z2 = self.nodes[l[2]][2]
        x3 = self.nodes[l[3]][0]
        y3 = self.nodes[l[3]][1]
        z3 = self.nodes[l[3]][2]
        # face 0
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
        # face 1
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
        # face 2
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
        # face 3
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
      elif len(l) == 5:
        # Pyramid
        x0 = self.nodes[l[0]][0]
        y0 = self.nodes[l[0]][1]
        z0 = self.nodes[l[0]][2]
        x1 = self.nodes[l[1]][0]
        y1 = self.nodes[l[1]][1]
        z1 = self.nodes[l[1]][2]
        x2 = self.nodes[l[2]][0]
        y2 = self.nodes[l[2]][1]
        z2 = self.nodes[l[2]][2]
        x3 = self.nodes[l[3]][0]
        y3 = self.nodes[l[3]][1]
        z3 = self.nodes[l[3]][2]
        x4 = self.nodes[l[4]][0]
        y4 = self.nodes[l[4]][1]
        z4 = self.nodes[l[4]][2]
        # face 0
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
        # face 1
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
        # face 2
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
        # face 3
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
        # face 4
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
      elif len(l) == 6:
        # Prism
        x0 = self.nodes[l[0]][0]
        y0 = self.nodes[l[0]][1]
        z0 = self.nodes[l[0]][2]
        x1 = self.nodes[l[1]][0]
        y1 = self.nodes[l[1]][1]
        z1 = self.nodes[l[1]][2]
        x2 = self.nodes[l[2]][0]
        y2 = self.nodes[l[2]][1]
        z2 = self.nodes[l[2]][2]
        x3 = self.nodes[l[3]][0]
        y3 = self.nodes[l[3]][1]
        z3 = self.nodes[l[3]][2]
        x4 = self.nodes[l[4]][0]
        y4 = self.nodes[l[4]][1]
        z4 = self.nodes[l[4]][2]
        x5 = self.nodes[l[5]][0]
        y5 = self.nodes[l[5]][1]
        z5 = self.nodes[l[5]][2]
        # face 0
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
        # face 1
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x5,y5,z5))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
        # face 2
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x5,y5,z5))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
        # face 3
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
        # face 4
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x5,y5,z5))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 5
      elif len(l) == 8:
        # Prism
        x0 = self.nodes[l[0]][0]
        y0 = self.nodes[l[0]][1]
        z0 = self.nodes[l[0]][2]
        x1 = self.nodes[l[1]][0]
        y1 = self.nodes[l[1]][1]
        z1 = self.nodes[l[1]][2]
        x2 = self.nodes[l[2]][0]
        y2 = self.nodes[l[2]][1]
        z2 = self.nodes[l[2]][2]
        x3 = self.nodes[l[3]][0]
        y3 = self.nodes[l[3]][1]
        z3 = self.nodes[l[3]][2]
        x4 = self.nodes[l[4]][0]
        y4 = self.nodes[l[4]][1]
        z4 = self.nodes[l[4]][2]
        x5 = self.nodes[l[5]][0]
        y5 = self.nodes[l[5]][1]
        z5 = self.nodes[l[5]][2]
        x6 = self.nodes[l[6]][0]
        y6 = self.nodes[l[6]][1]
        z6 = self.nodes[l[6]][2]
        x7 = self.nodes[l[7]][0]
        y7 = self.nodes[l[7]][1]
        z7 = self.nodes[l[7]][2]
        # face 0
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x5,y5,z5))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
        # face 1
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x5,y5,z5))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x6,y6,z6))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
        # face 2
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x6,y6,z6))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x7,y7,z7))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
        # face 3
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x7,y7,z7))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
        # face 4
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x0,y0,z0))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x1,y1,z1))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x2,y2,z2))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x3,y3,z3))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
        # face 5
        file.write("gp = geom.poly.create(gset)\n")
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x7,y7,z7))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x6,y6,z6))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x5,y5,z5))
        file.write("geom.poly.add.node(gset,gp,vector({},{},{}))\n".format(x4,y4,z4))
        file.write("geom.poly.close(gset,gp)\n")
        nlines += 6
      else:
        print "incorrect shape !"
    
    file.write("end\n")
    file.write("@create_mesh")
    nlines += 2
    file.close()
    
    
    itasca.command("""
    define read_tmp
      local a = array.create({nl})
      file.open("~tmp.txt",0,1) 
      file.read(a,{nl}) 
      file.close()
      array.command(a)
    end
    @read_tmp
    """.format(nl=nlines))
    os.remove("~tmp.txt")
