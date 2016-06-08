/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
  \\/     M anipulation  |
  -------------------------------------------------------------------------------
  License
  This file is not part of OpenFOAM.

  OpenFOAM is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

  Application
  demIcoFoam
  \*---------------------------------------------------------------------------*/

#include "demIcoFoam.H"
#include <stdexcept>
demIcoFoam::demIcoFoam() {
  Info << 0 << endl;
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
  piso_ = new pisoControl(*mesh_);
  Info << mesh_->nCells() << endl;

  // does it matter if this goes out of scope?
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
  f_ = new volVectorField(IOobject ("f", runTime_->timeName(), *mesh_,
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
  gradp_ = new volVectorField(fvc:: grad(*p_));
  cumulativeContErr_ = 0.0;
}

demIcoFoam::~demIcoFoam() {
  Info << "cleaning up demIcoFoam" << nl << endl;
  if (gradp_) delete gradp_;
  if (phi_) delete phi_;
  if (n_) delete n_;
  if (f_) delete f_;
  if (U_) delete U_;
  if (p_) delete p_;
  if (rho_) delete rho_;
  if (nu_) delete nu_;
  if (piso_) delete piso_;
  if (mesh_) delete mesh_;
  if (runTime_) delete runTime_;
  if (args_) delete args_;

}

void demIcoFoam::run(double time_increment) {
  volVectorField        &U = *U_;
  volVectorField        &f = *f_;
  volScalarField        &n = *n_;
  volScalarField        &p = *p_;
  dimensionedScalar     &nu = *nu_;
  surfaceScalarField    &phi = *phi_;
  Foam::Time            &runTime = *runTime_;
  Foam::fvMesh          &mesh = *mesh_;
  scalar &cumulativeContErr = cumulativeContErr_;

  runTime.setEndTime(runTime.value() + time_increment);
  while (runTime.loop())
  {
    Info<< "Time = " << runTime.timeName() << nl << endl;
    #include "CourantNo.H"

    fvVectorMatrix UEqn
      (fvm::ddt(U) + fvm::div(phi, U) -
       fvm::laplacian(nu, U) - f/n);

    if (piso_->momentumPredictor())
      solve(UEqn == -fvc::grad(p));

    while (piso_->correct())
    {
      volScalarField rAU(1.0/UEqn.A());

      volVectorField HbyA("HbyA", U);
      HbyA = rAU*UEqn.H();
      surfaceScalarField phiHbyA
        ("phiHbyA", (fvc::interpolate(HbyA) & mesh_->Sf()));

      adjustPhi(phiHbyA, U, p);

      while (piso_->correctNonOrthogonal())
      {
        fvScalarMatrix pEqn
          (fvm::laplacian(rAU, p) ==
           fvc::div(phiHbyA) //+ dndt
            );
        pEqn.setReference(pRefCell_, pRefValue_);
        pEqn.solve(mesh_->solver(p.select(piso_->finalInnerIter())));
        if (piso_->finalNonOrthogonalIter()) phi = phiHbyA - pEqn.flux();
      }


#include "continuityErrs.H"

      U = HbyA - rAU*fvc::grad(p);
      U.correctBoundaryConditions();
    }
    runTime_->write();
    Info<< "ExecutionTime = " << runTime_->elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime_->elapsedClockTime() << " s"
        << nl << endl;
  }
}


double demIcoFoam::flux_on_patch(char *patch_name)
{
  label inletPatchi = (*mesh_).boundaryMesh().findPatchID(patch_name);
  if (inletPatchi == -1)
    throw std::runtime_error("Cannot find boundary patch");
  scalar massFlux = sum((*phi_).boundaryField()[inletPatchi]);

  return massFlux;

}



//   volScalarField dndt = (n-oldn)/runTime.deltaT();


//   writeVectorField("u.dat", U);
//   writeScalarField("p.dat", p);
//   volVectorField gradp = fvc::grad(p);
//   writeVectorField("gradp.dat", gradp);
