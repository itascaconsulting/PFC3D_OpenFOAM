# distutils: language = c++
# distutils: sources = demIcoFoam.C
import numpy as np
cdef extern from "demIcoFoam.H":
   cdef cppclass demIcoFoam:
       demIcoFoam() except +
       int nCells()
       double rho()
       double nu()
       int nNodes()
       int nFaces()
       int face_node(int face, int node)
       int cell_face(int cell, int face)
       double face_center(int face, int j)
       double node_pos(int i, int j)
       int element(int i, int j)
       double f(int i, int j)
       void set_f(int i, int j, double v)
       double n(int i)
       void set_n(int i, double v)
       double p(int i)
       double U(int i, int j)
       double gradp(int i, int j)
       double phi(int face)
       void set_dt(double v)
       void run(double t) except +

cdef class pyDemIcoFoam:
    cdef demIcoFoam *thisptr
    def __cinit__(self): self.thisptr = new demIcoFoam()
    def __dealloc__(self): del self.thisptr
    def nCells(self): return self.thisptr.nCells()
    def rho(self): return self.thisptr.rho()
    def nu(self): return self.thisptr.nu()
    def mu(self): return self.nu()*self.rho()
    def nNodes(self): return self.thisptr.nNodes()
    def nFaces(self): return self.thisptr.nFaces()
    def set_dt(self, v): self.thisptr.set_dt(v)

    def faces(self):
        return np.array([[self.thisptr.face_node(i,j) for j in range(4)]
                         for i in range(self.nFaces())])

    def cell_faces(self):
        return np.array([[self.thisptr.cell_face(i,j) for j in range(6)]
                         for i in range(self.nCells())])

    def face_centers(self):
        return np.array([[self.thisptr.face_center(i,j) for j in range(3)]
                         for i in range(self.nFaces())])


    def nodes(self):
        return np.array([[self.thisptr.node_pos(i,0),
                          self.thisptr.node_pos(i,1),
                          self.thisptr.node_pos(i,2)]
                         for i in range(self.nNodes())])
    def elements(self):
        return np.array([[self.thisptr.element(i,j) for j in range(8)]
                         for i in range(self.nCells())])

    def f(self, value=None):
        if value is None:
            return np.array([[self.thisptr.f(i,0),
                              self.thisptr.f(i,1),
                              self.thisptr.f(i,2)]
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(), 3)
            for i in range(self.nCells()):
                for j in range(3):
                    self.thisptr.set_f(i,j,value[i][j])

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
