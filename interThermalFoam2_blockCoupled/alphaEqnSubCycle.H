label nAlphaCorr
(
    readLabel(pimple.lookup("nAlphaCorr"))
);

label nAlphaSubCycles
(
    readLabel(pimple.lookup("nAlphaSubCycles"))
);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum = 0.0*rhoPhi;
    
    surfaceScalarField rhocpPhi1Sum = 0.0*rhocpPhi1;
    surfaceScalarField rhocpPhi2Sum = 0.0*rhocpPhi2;

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
#       include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;

        rhocpPhi1Sum += (runTime.deltaT()/totalDeltaT)*rhocpPhi1;
        rhocpPhi2Sum += (runTime.deltaT()/totalDeltaT)*rhocpPhi2;

    }

    rhoPhi = rhoPhiSum;
    
    rhocpPhi1 = rhocpPhi1Sum;
    rhocpPhi2 = rhocpPhi2Sum;
}
else
{
#       include "alphaEqn.H"
}

interface.correct();

rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;

rhocpalpha1 == alpha1*rho1*cp1;
rhocpalpha2 == (1.0-alpha1)*rho2*cp2;
