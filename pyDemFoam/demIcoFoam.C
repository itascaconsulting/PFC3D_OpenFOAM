/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
  \\/     M anipulation  |
  ----------------------------------------------------------------------------
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
no
  Application
  demIcoFoam
  \*------------------------------------------------------------------------*/

#include "demIcoFoam.H"

demIcoFoam::demIcoFoam() {
  Info << "creating demIcoFoam object" << nl << endl;
  f_ = new volVectorField(IOobject ("f", runTime_->timeName(), *mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE), *mesh_);
  piso_ = new pisoControl(*mesh_);
}

demIcoFoam::~demIcoFoam() {
  Info << "cleaning up demIcoFoam" << nl << endl;
  if (piso_) delete piso_;
  if (f_) delete f_;
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
        pEqn.solve(mesh.solver(p.select(piso_->finalInnerIter())));
        if (piso_->finalNonOrthogonalIter()) phi = phiHbyA - pEqn.flux();
      }


#include "continuityErrs.H"

      U = HbyA - rAU*fvc::grad(p);
      U.correctBoundaryConditions();
    }
    runTime.write();
    Info<< "ExecutionTime = " << runTime_->elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime_->elapsedClockTime() << " s"
        << nl << endl;
  }
}
