/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
/*
    //Heat transfer relaxation:
    //Interphase transport model
    //Interface area density (m^2/m^3)
    const volVectorField gradAlpha = fvc::grad(alpha1);
    const volScalarField Aint = mag(gradAlpha);
    //Get normal to interface
    const dimensionedScalar deltaN("deltaN", dimensionSet(0,-1,0,0,0,0,0), 1E-16);
	const volVectorField ns = gradAlpha/(mag(gradAlpha) + deltaN);
	//Average k in each cell
	const volScalarField kc = (alpha1*k1 + alpha2*k2);
    //hA_effective in each cell
    const volScalarField hA_eff = kc*mag(gradAlpha)*Aint;
    
    const dimensionedScalar deltat("dt", dimensionSet(0,0,1,0,0,0,0), runTime.deltaTValue() );
    const volScalarField DeltaT1 = deltat*hA_eff*(T2-T1)/(rhocpalpha1);
    const volScalarField DeltaT2 = deltat*hA_eff*(T1-T2)/(rhocpalpha2);    
    const volScalarField DeltaTRatio = max( DeltaT1/(T2-T1), DeltaT2/(T1-T2) );
    const scalar RelaxationFactor = 0.4 / gMax( DeltaTRatio );
    Info<< endl << "Required relaxation factor: " << RelaxationFactor << endl << endl;
*/


    scalar maxDeltaTFact =
//        min( min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL)), RelaxationFactor );
        min(maxCo/(CoNum + SMALL), maxAlphaCo/(alphaCoNum + SMALL));

    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
