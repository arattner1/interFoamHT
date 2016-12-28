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

    surfaceScalarField rhocpPhiSum1
    (
        IOobject
        (
            "rhocpPhiSum1",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhocpPhi1.dimensions(), 0)
    );

    surfaceScalarField rhocpPhiSum2
    (
        IOobject
        (
            "rhocpPhiSum2",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhocpPhi1.dimensions(), 0)
    );


    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        rhocpPhiSum1 += (runTime.deltaT()/totalDeltaT)*rhocpPhi1;
        rhocpPhiSum2 += (runTime.deltaT()/totalDeltaT)*rhocpPhi2;
    }

    rhoPhi = rhoPhiSum;
    rhocpPhi1 = rhocpPhiSum1;
    rhocpPhi2 = rhocpPhiSum2;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;
rhocpalpha1 == alpha1*rho1*cp1;
rhocpalpha2 == alpha2*rho2*cp2;
