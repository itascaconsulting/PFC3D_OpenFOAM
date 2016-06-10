# distutils: language = c++
# distutils: sources = demIcoFoam.C

import numpy as np

from demBaseFoam cimport demBaseFoam

cdef class pyDemBaseFoam:
    cdef demBaseFoam *thisptr
    def __cinit__(self):
        # we should not allocate here?
        pass # self.thisptr = new demIcoFoam()

    def __dealloc__(self):
        # we do not want to free here?
        pass #del self.thisptr

    def nCells(self): return self.thisptr.nCells()
    def rho(self): return self.thisptr.rho()
    def nu(self): return self.thisptr.nu()
    def mu(self): return self.nu()*self.rho()
    def nNodes(self): return self.thisptr.nNodes()
    def nFaces(self): return self.thisptr.nFaces()
    def set_dt(self, v): self.thisptr.set_dt(v)
    def dt(self): return self.thisptr.dt()
    def cell_near(self, x,y,z): return self.thisptr.cell_near(x,y,z)
    def flux_on_patch(self, name): return self.thisptr.flux_on_patch(name)

    def faces(self):
        return np.array([[self.thisptr.face_node(i,j) for j in range(4)]
                         for i in range(self.nFaces())])

    def cell_faces(self):
        return np.array([[self.thisptr.cell_face(i,j) for j in range(6)]
                         for i in range(self.nCells())])

    def face_centers(self):
        return np.array([[self.thisptr.face_center(i,j) for j in range(3)]
                         for i in range(self.nFaces())])

    def cell_centers(self):
        return np.array([[self.thisptr.cell_center(i,j) for j in range(3)]
                         for i in range(self.nCells())])

    def cell_volumes(self):
        return np.array([self.thisptr.cell_volume(i)
                         for i in range(self.nCells())])

    def cell_fluxes(self):
        return np.array([[self.thisptr.cell_flux(i,j) for j in range(6)]
                         for i in range(self.nCells())])

    def nodes(self):
        return np.array([[self.thisptr.node_pos(i,0),
                          self.thisptr.node_pos(i,1),
                          self.thisptr.node_pos(i,2)]
                         for i in range(self.nNodes())])
    def elements(self):
        return np.array([[self.thisptr.element(i,j) for j in range(8)]
                         for i in range(self.nCells())])

    def n(self, value=None):
        if value is None:
            return np.array([self.thisptr.n(i)
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(),)
            for i in range(self.nCells()):
                self.thisptr.set_n(i,value[i])

    def p(self):
        return np.array([self.thisptr.p(i)
                         for i in range(self.nCells())])

    def phi(self):
        return np.array([self.thisptr.phi(i)
                         for i in range(self.nFaces())])

    def U(self):
        return np.array([[self.thisptr.U(i,0),
                          self.thisptr.U(i,1),
                          self.thisptr.U(i,2)]
                         for i in range(self.nCells())])

    def gradp(self):
        return np.array([[self.thisptr.gradp(i,0),
                          self.thisptr.gradp(i,1),
                          self.thisptr.gradp(i,2)]
                         for i in range(self.nCells())])

    def solve(self, time_increment):
        self.thisptr.run(time_increment)


cdef extern from "demIcoFoam.H":
   cdef cppclass demIcoFoam(demBaseFoam):
       demIcoFoam() except +
       double f(int i, int j)
       double set_f(int i, int j, double value)


cdef class pyDemIcoFoam(pyDemBaseFoam):
    cdef demIcoFoam *dthisptr
    def __cinit__(self):
        self.dthisptr = self.thisptr = new demIcoFoam()
    def __dealloc__(self):
        del self.thisptr

    def f(self, value=None):
        if value is None:
            return np.array([[self.dthisptr.f(i,0),
                              self.dthisptr.f(i,1),
                              self.dthisptr.f(i,2)]
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(), 3)
            for i in range(self.nCells()):
                for j in range(3):
                    self.dthisptr.set_f(i,j,value[i][j])
