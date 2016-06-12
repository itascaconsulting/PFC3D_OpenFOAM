# distutils: language = c++
# distutils: sources = demIcoFoam.C demSimpleFoam.C demBaseFoam.C

import numpy as np

cdef extern from "demBaseFoam.H":
   cdef cppclass demBaseFoam:
       demBaseFoam() except +
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
       double n(int i)
       void set_n(int i, double v)
       double p(int i)
       double U(int i, int j)
       double gradp(int i, int j)
       double phi(int face) except +
       void set_dt(double v)
       double dt()
       int cell_near(double x, double y, double z)
       double cell_center(int cell, int j)
       double cell_volume(int cell)
       double flux_on_patch(char *patch_name) except +
       void run(double t) except +


cdef class pyDemBaseFoam:
    cdef demBaseFoam *thisptr
    def __cinit__(self): pass
    def __dealloc__(self): pass
    def nCells(self):
        """() -> int. Return the number of cells in the CFD mesh. """
        return self.thisptr.nCells()
    def rho(self):
        """() -> float. Return the fluid density."""
        return self.thisptr.rho()
    def nu(self):
        """() -> float. Return the fluid kinematic viscosity."""
        return self.thisptr.nu()
    def mu(self):
        """() -> float. Return the fluid dynamic viscosity."""
        return self.nu()*self.rho()
    def nNodes(self):
        """() -> int. Return the number of nodes in the CFD mesh."""
        return self.thisptr.nNodes()
    def nFaces(self):
        """() -> int. Return the number of faces in the CFD mesh."""
        return self.thisptr.nFaces()
    def set_dt(self, v):
        """(time_step: float) -> None. Set the OpenFOAM timestep."""
        self.thisptr.set_dt(v)
    def dt(self):
        """() -> float. Return the OpenFOAM timestep."""
        return self.thisptr.dt()
    def cell_near(self, x,y,z):
        """(x: float, y: float, z: float) -> int. Return the CFD mesh cell
nearst the given point. Returns -1 if a call cannot be found."""
        return self.thisptr.cell_near(x,y,z)
    def flux_on_patch(self, name):
        """(patch: string) -> float. Return the mass flux on the given boundary patch."""
        return self.thisptr.flux_on_patch(name)

    def faces(self):
        """() -> int array{n, 4}. Return an array of integers which describe
the faces of the elements. n is the number of faces. The integers are
indices into the nodes array.
        """
        return np.array([[self.thisptr.face_node(i,j) for j in range(4)]
                         for i in range(self.nFaces())])

    def cell_faces(self):
        return np.array([[self.thisptr.cell_face(i,j) for j in range(6)]
                         for i in range(self.nCells())])

    def face_centers(self):
        """() -> float array{n,3}. Return an array giving the center of each face. n is the number of faces."""
        return np.array([[self.thisptr.face_center(i,j) for j in range(3)]
                         for i in range(self.nFaces())])

    def cell_centers(self):
        """() -> float array{n,3}. Return an array of the center of each cell. n is the number of cells."""
        return np.array([[self.thisptr.cell_center(i,j) for j in range(3)]
                         for i in range(self.nCells())])

    def cell_volumes(self):
        """() -> float array{n}. Return the volume of the cells."""
        return np.array([self.thisptr.cell_volume(i)
                         for i in range(self.nCells())])

    def cell_fluxes(self):
        return np.array([[self.thisptr.cell_flux(i,j) for j in range(6)]
                         for i in range(self.nCells())])

    def nodes(self):
        """() -> array{n,3}. Return a numpy array of the CFD mesh node locations where n is the number of nodes in the mesh."""
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



########################################################
## derived
########################################################

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

###########################

cdef extern from "demSimpleFoam.H":
   cdef cppclass demSimpleFoam(demBaseFoam):
       demSimpleFoam() except +
       double ubar(int i, int j)
       void set_ubar(int i, int j, double v)
       double beta(int i)
       void set_beta(int i, double v)


cdef class pyDemSimpleFoam(pyDemBaseFoam):
    cdef demSimpleFoam *dthisptr
    def __cinit__(self):
        self.dthisptr = self.thisptr = new demSimpleFoam()
    def __dealloc__(self):
        del self.thisptr

    def ubar(self, value=None):
        if value is None:
            return np.array([[self.dthisptr.ubar(i,0),
                              self.dthisptr.ubar(i,1),
                              self.dthisptr.ubar(i,2)]
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(), 3)
            for i in range(self.nCells()):
                for j in range(3):
                    self.dthisptr.set_ubar(i,j,value[i][j])

    def beta(self, value=None):
        if value is None:
            return np.array([self.dthisptr.beta(i)
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(),)
            for i in range(self.nCells()):
                self.dthisptr.set_beta(i,value[i])

cdef extern from "demSimpleFoam2.H":
   cdef cppclass demSimpleFoam2(demBaseFoam):
       demSimpleFoam2() except +
       double ubar(int i, int j)
       void set_ubar(int i, int j, double v)
       double beta(int i)
       void set_beta(int i, double v)


cdef class pyDemSimpleFoam2(pyDemBaseFoam):
    cdef demSimpleFoam2 *dthisptr
    def __cinit__(self):
        self.dthisptr = self.thisptr = new demSimpleFoam2()
    def __dealloc__(self):
        del self.thisptr

    def ubar(self, value=None):
        if value is None:
            return np.array([[self.dthisptr.ubar(i,0),
                              self.dthisptr.ubar(i,1),
                              self.dthisptr.ubar(i,2)]
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(), 3)
            for i in range(self.nCells()):
                for j in range(3):
                    self.dthisptr.set_ubar(i,j,value[i][j])

    def beta(self, value=None):
        if value is None:
            return np.array([self.dthisptr.beta(i)
                             for i in range(self.nCells())])
        else:
            value = np.asarray(value, dtype=np.double)
            assert value.shape == (self.nCells(),)
            for i in range(self.nCells()):
                self.dthisptr.set_beta(i,value[i])
