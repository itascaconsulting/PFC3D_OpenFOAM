/*---------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
  \\/     M anipulation  |
  ---------------------------------------------------------------------
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
  demSimpleFoam
  \*--------------------------------------------------------------*/

#include "demSimpleFoam.H"

demSimpleFoam::demSimpleFoam() {
  Info << "creating demSimpleFoam object" << nl << endl;

  ubar_ = new volVectorField(IOobject ("ubar", runTime_->timeName(), *mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE), *mesh_);
  beta_ = new volScalarField (IOobject ("beta", runTime_->timeName(), *mesh_,
                                     IOobject::MUST_READ,
                                     IOobject::AUTO_WRITE), *mesh_);
  simple_ = new simpleControl(*mesh_);
}

demSimpleFoam::~demSimpleFoam() {
  Info << "cleaning up demSimpleFoam" << nl << endl;
  if (ubar_) delete ubar_;
  if (beta_) delete beta_;
}

void demSimpleFoam::run(double) {
  volScalarField &beta = *beta_;
  volVectorField &ubar = *ubar_;
  volVectorField &U = *U_;
  volScalarField &n = *n_;
  volScalarField        &p = *p_;
  dimensionedScalar     &nu = *nu_;
  surfaceScalarField    &phi = *phi_;
  Foam::Time            &runTime = *runTime_;
  Foam::fvMesh          &mesh = *mesh_;
  scalar &cumulativeContErr = cumulativeContErr_;

  while (simple_->loop())
  {
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // --- Pressure-velocity SIMPLE corrector
    {
      // Momentum predictor
      tmp<fvVectorMatrix> UEqn (fvm::div(phi, U) -
                                fvm::laplacian(nu, U) -
                                beta*ubar/n -
                                fvm::Sp(beta/n, U)
        );

      UEqn().relax();

      solve(UEqn() == -fvc::grad(p));

      //#include "pEqn.H"
      {
        volScalarField rAU(1.0/UEqn().A());
        volVectorField HbyA("HbyA", U);
        HbyA = rAU*UEqn().H();

        surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA) &
                                   mesh.Sf());
        adjustPhi(phiHbyA, U, p);

        tmp<volScalarField> rAtU(rAU);

        if (simple_->consistent())
        {
          rAtU = 1.0/(1.0/rAU - UEqn().H1());
          phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
          HbyA -= (rAU - rAtU())*fvc::grad(p);
        }

        UEqn.clear();

        // Non-orthogonal pressure corrector loop
        while (simple_->correctNonOrthogonal())
        {
          fvScalarMatrix pEqn (fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA));

          pEqn.setReference(pRefCell_, pRefValue_);

          pEqn.solve();

          if (simple_->finalNonOrthogonalIter())
          {
            phi = phiHbyA - pEqn.flux();
          }
        }

        #include "continuityErrs.H"

        // Explicitly relax pressure for momentum corrector
        p.relax();

        // Momentum corrector
        U = HbyA - rAtU()*fvc::grad(p);
        U.correctBoundaryConditions();
      }
    }
    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
  }
  Info << "SIMPLE Solve Ended. \n" << endl;
}
