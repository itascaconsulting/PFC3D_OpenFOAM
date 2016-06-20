#pragma once
#include "fvCFD.H"
#include "pisoControl.H"
#include "meshSearch.H"

#include <iostream>
#include <fstream>
#include <stdexcept>

#define VEC_GET(T) double T(int i, int j) { return (*T##_)[i][j]; }
#define VEC_SET(T) void set_##T(int i, int j, double v) { (*T##_)[i][j]=v; }
#define VEC_SETGET(T) VEC_SET(T) VEC_GET(T)

#define GET(T) double T(int i) { return (*T##_)[i]; }
#define SET(T) void set_##T(int i, double v) { (*T##_)[i]=v; }
#define SETGET(T) SET(T) GET(T)

class demBaseFoam {
public:
  demBaseFoam();
  virtual ~demBaseFoam();
  int nCells() { return mesh_->nCells(); }
  double rho() { return rho_->value(); }
  double nu() { return nu_->value(); }
  int nNodes() { return mesh_->nPoints(); }
  int nFaces() { return mesh_->nFaces(); }
  int face_node(int face, int node) { return mesh_->faces()[face][node]; }
  int cell_face(int cell, int face) { return mesh_->cells()[cell][face]; }
  double cell_flux(int cell, int face);
  int cell_near(double x, double y, double z);
  double flux_on_patch(char *patch_name);
  double cell_center(int cell, int j) { return mesh_->cellCentres()[cell][j]; }
  double cell_volume(int cell) { return mesh_->cellVolumes()[cell]; }
  double face_center(int face, int j) { return mesh_->faceCentres()[face][j]; }
  double node_pos(int i, int j) { return mesh_->points()[i][j]; }
  int element(int i, int j) { return mesh_->cellShapes()[i][j]; }
  void set_dt(double t) { runTime_->setDeltaT(t); }
  double dt() { return runTime_->deltaT().value(); }
  void set_endtime(double t) { runTime_->setEndTime(t); }
  virtual void run(double v)=0;

  SETGET(n)
  VEC_SETGET(U)
  SETGET(p)
  VEC_GET(gradp)
  GET(phi)

  Foam::argList *args_;
  Foam::Time  *runTime_;
  Foam::fvMesh *mesh_;

  dimensionedScalar *nu_;
  dimensionedScalar *rho_;
  volScalarField *p_;
  volVectorField *U_;
  volScalarField *n_;
  surfaceScalarField *phi_;
  label pRefCell_;
  scalar pRefValue_;

  volVectorField *gradp_;
  scalar cumulativeContErr_;
};
