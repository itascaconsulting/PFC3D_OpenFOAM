# distutils: language = c++
# distutils: sources = demBaseFoam.C
# -*- cython-mode -*-

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
