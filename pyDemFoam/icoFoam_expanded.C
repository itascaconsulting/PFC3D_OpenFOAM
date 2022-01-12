
#include "fvCFD.H"
#include "pisoControl.H"


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

//#include "postProcess.H" might not need this as this mostly reads
// command line argument and looks at the results of previous models.

//#include "addCheckCaseOptions.H" // looks like we do not need this one
// #include "setRootCaseLists.H" // or this on these are related to command line args

    // #include "createTime.H"
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time runTime(Foam::Time::controlDictName, args);

    //#include "createMesh.H"
    Foam::fvMesh mesh (Foam::IOobject(Foam::fvMesh::defaultRegion, runTime.timeName(), runTime, Foam::IOobject::MUST_READ));

    pisoControl piso(mesh);

    // #include "createFields.H"
    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties(IOobject("transportProperties", runTime.constant(), mesh, IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE));

    dimensionedScalar nu ("nu", dimViscosity, transportProperties);

    Info<< "Reading field p\n" << endl;
    volScalarField p (IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);


    Info<< "Reading field U\n" << endl;
    volVectorField U (IOobject ("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);


   //#include "createPhi.H"
    Info<< "Reading/calculating face flux field phi\n" << endl;

    surfaceScalarField phi (IOobject ("phi", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), fvc::flux(U));



    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());


//#include "initContinuityErrs.H"
    //uniformDimensionedScalarField cumulativeContErrIO (IOobject ("cumulativeContErr", runTime.timeName(), "uniform", mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), dimensionedScalar(dimless, Zero));
    // I think this is just supposed to init to zero in the absence of other run data
    scalar& cumulativeContErr = 0;cumulativeContErrIO.value();


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#include "CourantNo.H" // not sure which actual file is included here
// I think it is this?? OpenFOAM-v2112/src/finiteVolume/cfdTools/incompressible/CourantNo.H
        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

#include "continuityErrs.H"
// OpenFOAM-v2112/src/finiteVolume/cfdTools/incompressible/continuityErrs.H
// I think this is the file, should be OK to leave as an include
            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
