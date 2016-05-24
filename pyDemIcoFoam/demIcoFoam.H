#include "fvCFD.H"
#include "pisoControl.H"
#include <iostream>
#include <fstream>

#define VEC_GET(T) double T(int i, int j) { return (*T##_)[i][j]; }
#define VEC_SET(T) void set_##T(int i, int j, double v) { (*T##_)[i][j]=v; }
#define VEC_SETGET(T) VEC_SET(T) VEC_GET(T)

#define GET(T) double T(int i) { return (*T##_)[i]; }
#define SET(T) void set_##T(int i, double v) { (*T##_)[i]=v; }
#define SETGET(T) SET(T) GET(T)


class demIcoFoam {
public:
  demIcoFoam();
  virtual ~demIcoFoam() { Info << "cleaning up" << nl << endl;}
  int nCells() { return mesh_->nCells(); }
  double rho() { return rho_->value(); }
  double nu() { return nu_->value(); }
  int nNodes() { return mesh_->nPoints(); }
  double node_pos(int i, int j) { return mesh_->points()[i][j]; }
  int element(int i, int j) { return mesh_->cellShapes()[i][j]; }
  void set_dt(double t) { runTime_->setDeltaT(t); }
  void set_endtime(double t) { runTime_->setEndTime(t); }
  void run(double v);
  VEC_SETGET(f)
  SETGET(n)
  VEC_SETGET(U)
  SETGET(p)
  VEC_GET(gradp)

  Foam::argList *args_;
  Foam::Time  *runTime_;
  Foam::fvMesh *mesh_;
  pisoControl  *piso_;
  dimensionedScalar *nu_;
  dimensionedScalar *rho_;
  volScalarField *p_;
  volVectorField *U_;
  volVectorField *f_;
  volScalarField *n_;
  surfaceScalarField *phi_;
  label pRefCell_;
  scalar pRefValue_;

  volVectorField *gradp_;
  scalar cumulativeContErr_;
};
