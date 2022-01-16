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

#include "demIcoFoamSemiImplicitDrag.H"

demIcoFoamSemiImplicitDrag::demIcoFoamSemiImplicitDrag() {
  Info << "creating demIcoFoamSemiImplicitDrag object" << nl << endl;
  ubar_ = new volVectorField(IOobject ("ubar", runTime_->timeName(), *mesh_,
                                    IOobject::MUST_READ,
                                    IOobject::AUTO_WRITE), *mesh_);
  beta_ = new volScalarField (IOobject ("beta", runTime_->timeName(), *mesh_,
                                     IOobject::MUST_READ,
                                     IOobject::AUTO_WRITE), *mesh_);
  piso_ = new pisoControl(*mesh_);
}

demIcoFoamSemiImplicitDrag::~demIcoFoamSemiImplicitDrag() {
  Info << "cleaning up demIcoFoamSemiImplicitDrag" << nl << endl;
  if (piso_) delete piso_;
  if (ubar_) delete ubar_;
  if (beta_) delete beta_;
}

void demIcoFoamSemiImplicitDrag::run(double time_increment) {
  volScalarField &beta = *beta_; // aliases for convenience
  volVectorField &ubar = *ubar_;
  volVectorField        &U = *U_;
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
      (fvm::ddt(U)
       + fvm::div(phi, U)
       - fvm::laplacian(nu, U)
       - beta*ubar/n
       + fvm::Sp(beta/n, U));

    if (piso_->momentumPredictor())
      solve(UEqn == -fvc::grad(p));

    while (piso_->correct())
    {
      volScalarField rAU(1.0/UEqn.A());
      // volVectorField HbyA("HbyA", U);
      // HbyA = rAU*UEqn.H();
      // surfaceScalarField phiHbyA
      //   ("phiHbyA", (fvc::interpolate(HbyA) & mesh_->Sf()));

      volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
      surfaceScalarField phiHbyA ("phiHbyA", fvc::flux(HbyA) + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi));

      adjustPhi(phiHbyA, U, p);

      // Update the pressure BCs to ensure flux consistency
      constrainPressure(p, U, phiHbyA, rAU);


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
