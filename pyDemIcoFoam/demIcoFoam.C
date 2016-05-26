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
  runTime_->setEndTime(runTime_->value() + time_increment);
  while (runTime_->loop())
  {
    {
      Info<< "Time = " << runTime_->timeName() << nl << endl;
      scalar CoNum = 0.0;
      scalar meanCoNum = 0.0;
      if (mesh_->nInternalFaces())
      {
        scalarField sumPhi (fvc::surfaceSum(mag(*phi_))().internalField());
        CoNum = 0.5*gMax(sumPhi/mesh_->V().field())*runTime_->deltaTValue();
        meanCoNum =
          0.5*(gSum(sumPhi)/gSum(mesh_->V().field()))*runTime_->deltaTValue();
      }
      Info<< "Courant Number mean: " << meanCoNum
          << " max: " << CoNum << endl;
    }


    fvVectorMatrix UEqn
      (fvm::ddt(*n_,*U_) + (*n_)*fvm::div(*phi_, *U_) -
       fvm::laplacian(*nu_, *U_) - *f_);

    if (piso_->momentumPredictor())
      solve(UEqn == -(*n_)*fvc::grad(*p_));

    while (piso_->correct())
    {
      volScalarField rAU(1.0/UEqn.A());

      volVectorField HbyA("HbyA", *U_);
      HbyA = rAU*UEqn.H();
      surfaceScalarField phiHbyA
        ("phiHbyA", (fvc::interpolate(HbyA) & mesh_->Sf()));

      adjustPhi(phiHbyA, *U_, *p_);

      while (piso_->correctNonOrthogonal())
      {
        fvScalarMatrix pEqn
          (fvm::laplacian(rAU*(*n_)*(*n_), *p_) ==
           fvc::div(fvc::interpolate(*n_)*phiHbyA) //+ dndt
            );
        pEqn.setReference(pRefCell_, pRefValue_);
        pEqn.solve(mesh_->solver(p_->select(piso_->finalInnerIter())));
        if (piso_->finalNonOrthogonalIter()) *phi_ = phiHbyA - pEqn.flux();
      }

      {
        //volScalarField contErr(fvc::div(fvc::interpolate(n)*phi)+dndt);
        volScalarField contErr(fvc::div(fvc::interpolate(*n_)*(*phi_)));

        scalar sumLocalContErr = runTime_->deltaTValue()*
          mag(contErr)().weightedAverage(mesh_->V()).value();

        scalar globalContErr = runTime_->deltaTValue()*
          contErr.weightedAverage(mesh_->V()).value();
        cumulativeContErr_ += globalContErr;

        Info<< "time step continuity errors : sum local = " << sumLocalContErr
            << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr_
            << endl;
      }
      *U_ = HbyA - (*n_)*rAU*fvc::grad(*p_);
      U_->correctBoundaryConditions();
    }
    runTime_->write();
    Info<< "ExecutionTime = " << runTime_->elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime_->elapsedClockTime() << " s"
        << nl << endl;
  }
}





//   volScalarField dndt = (n-oldn)/runTime.deltaT();


//   writeVectorField("u.dat", U);
//   writeScalarField("p.dat", p);
//   volVectorField gradp = fvc::grad(p);
//   writeVectorField("gradp.dat", gradp);
