if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhoPhi.dimensions(), 0)
    );

    surfaceScalarField rhocpPhiSum
    (
        IOobject
        (
            "rhocpPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhocpPhi.dimensions(), 0)
    );



    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        rhocpPhiSum += (runTime.deltaT()/totalDeltaT)*rhocpPhi;
    }

    rhoPhi = rhoPhiSum;
    rhocpPhi = rhocpPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;
rhocp == alpha1*rho1*cp1 + alpha2*rho2*cp2;
