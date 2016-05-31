# distutils: language = c++
# distutils: sources = demSimpleFoam.C
import numpy as np
cdef extern from "demSimpleFoam.H":
   cdef cppclass demSimpleFoam:
       demSimpleFoam() except +
       int nCells()
       double rho()
       double nu()
       int nNodes()
       int nFaces()
       int face_node(int face, int node)
       int cell_face(int cell, int face)
       double face_center(int face, int j)
       double cell_flux(int cell, int face)
       double node_pos(int i, int j)
       int element(int i, int j)
       double ubar(int i, int j)
       void set_ubar(int i, int j, double v)
       double beta(int i)
       void set_beta(int i, double v)
       double n(int i)
       void set_n(int i, double v)
       double p(int i)
       double U(int i, int j)
       double gradp(int i, int j)
       double phi(int face) except +
       void run() except +
       int cell_near(double x, double y, double z)
       double cell_center(int cell, int j)
       double cell_volume(int cell)
       double flux_on_patch(char *patch_name) except +


cdef class pyDemSimpleFoam:
    cdef demSimpleFoam *thisptr
    def __cinit__(self): self.thisptr = new demSimpleFoam()
    def __dealloc__(self): del self.thisptr
    def nCells(self): return self.thisptr.nCells()
    def rho(self): return self.thisptr.rho()
    def nu(self): return self.thisptr.nu()
    def mu(self): return self.nu()*self.rho()
    def nNodes(self): return self.thisptr.nNodes()
    def nFaces(self): return self.thisptr.nFaces()
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

    def ubar(self, value=None):
        if value is None:
            return np.array([[self.thisptr.ubar(i,0),
                              self.thisptr.ubar(i,1),
                              self.thisptr.ubar(i,2)]
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(), 3)
            for i in range(self.nCells()):
                for j in range(3):
                    self.thisptr.set_ubar(i,j,value[i][j])

    def n(self, value=None):
        if value is None:
            return np.array([self.thisptr.n(i)
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(),)
            for i in range(self.nCells()):
                self.thisptr.set_n(i,value[i])

    def beta(self, value=None):
        if value is None:
            return np.array([self.thisptr.beta(i)
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(),)
            for i in range(self.nCells()):
                self.thisptr.set_beta(i,value[i])


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

    def solve(self):
        self.thisptr.run()
