#include "demBaseFoam.H"

demBaseFoam::demBaseFoam() {
  int argc=1;
  const char *argv[]={"demIcoFoam", NULL};
  char ** argv2 = const_cast<char **>(argv);
  args_ = new Foam::argList(argc, argv2);
  if (!args_->checkRootCase()) Foam::FatalError.exit();
  runTime_ = new Foam::Time(Foam::Time::controlDictName, *args_);
  mesh_ = new Foam::fvMesh(Foam::IOobject (Foam::fvMesh::defaultRegion,
                                           runTime_->timeName(),
                                           *runTime_,
                                           Foam::IOobject::MUST_READ));


  IOdictionary transportProperties(IOobject("transportProperties",
                                            runTime_->constant(),
                                            *mesh_,
                                            IOobject::MUST_READ_IF_MODIFIED,
                                            IOobject::NO_WRITE));

  nu_ = new dimensionedScalar("nu", dimViscosity,
                              transportProperties.lookup("nu"));
  rho_ = new dimensionedScalar("rho", dimDensity,
                               transportProperties.lookup("rho"));

  p_ = new volScalarField (IOobject ("p", runTime_->timeName(), *mesh_,
                                     IOobject::MUST_READ,
                                     IOobject::AUTO_WRITE), *mesh_);
  U_ = new volVectorField (IOobject ("U", runTime_->timeName(), *mesh_,
                                     IOobject::MUST_READ,
                                     IOobject::AUTO_WRITE), *mesh_);
  n_ = new volScalarField(IOobject ("n", runTime_->timeName(), *mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE), *mesh_);
  phi_ = new surfaceScalarField (IOobject ("phi", runTime_->timeName(), *mesh_,
                                           IOobject::READ_IF_PRESENT,
                                           IOobject::AUTO_WRITE),
                                 linearInterpolate(*U_) & mesh_->Sf());

  pRefCell_ = 0;
  pRefValue_ = 0.0;
  setRefCell(*p_, mesh_->solutionDict().subDict("PISO"), pRefCell_, pRefValue_);
  mesh_->setFluxRequired(p_->name());
  gradp_ = new volVectorField(fvc::grad(*p_));
  cumulativeContErr_ = 0.0;
}

demBaseFoam::~demBaseFoam() {
  Info << "cleaning up demBaseFoam" << endl;
  if (gradp_) delete gradp_;
  if (phi_) delete phi_;
  if (n_) delete n_;
  if (U_) delete U_;
  if (p_) delete p_;
  if (rho_) delete rho_;
  if (nu_) delete nu_;
  if (mesh_) delete mesh_;
  if (runTime_) delete runTime_;
  if (args_) delete args_;
}

double demBaseFoam::flux_on_patch(char *patch_name)
{
  label inletPatchi = (*mesh_).boundaryMesh().findPatchID(patch_name);
  if (inletPatchi == -1)
    throw std::runtime_error("Cannot find boundary patch");
  scalar massFlux = sum((*phi_).boundaryField()[inletPatchi]);
  return massFlux;
}

double demBaseFoam::cell_flux(int cell, int face) {
  label gface = mesh_->cells()[cell][face];
  if (mesh_->faceOwner()[gface]==cell) return (*phi_)[gface];
  else return -(*phi_)[gface];
}

int demBaseFoam::cell_near(double x, double y, double z) {
  meshSearch search(*mesh_);
  return search.findNearestCell(point(x,y,z));
}
