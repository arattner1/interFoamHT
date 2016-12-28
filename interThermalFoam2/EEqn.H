    //Limiting to keep problem conditioned, otherwise the equations get a 
    //divided by 0 type issue
    const scalar eps = 1E-8;
    rhocpalpha1 = max( rhocpalpha1, eps*rho1*cp1);
    rhocpalpha2 = max( rhocpalpha2, eps*rho2*cp2);


    //Get conductivity on cell faces
//    const surfaceScalarField alpha1f = fvc::interpolate(alpha1);
//    const surfaceScalarField kf1 = alpha1f*k1;
//    const surfaceScalarField kf2 = (1.0-alpha1f)*k2;


    //Interphase transport model
    //Interface area density (m^2/m^3)
    const volVectorField gradAlpha = fvc::grad(alpha1);
    const volScalarField Aint = mag(gradAlpha)*pos(alpha1-0.01)*pos(alpha2-0.01);
    //Get normal to interface
    const dimensionedScalar deltaN("deltaN", dimensionSet(0,-1,0,0,0,0,0), 1E-16);
	const volVectorField ns = gradAlpha/(mag(gradAlpha) + deltaN);
	//Average k in each cell
	const volScalarField kc = (alpha1*k1 + alpha2*k2);
    //Inter-phase heat transfer
    const volScalarField Q12 = -kc*(fvc::grad(T)&ns)*Aint;

    fvScalarMatrix EEqn1
    (
        fvm::ddt(rhocpalpha1, T1)
      + fvm::div(rhocpPhi1, T1)
//      - fvm::laplacian(kf1, T1)
      - fvm::laplacian(k1, T1)
      == Q12
    );

    EEqn1.solve();

    fvScalarMatrix EEqn2
    (
        fvm::ddt(rhocpalpha2, T2)
      + fvm::div(rhocpPhi2, T2)
//      - fvm::laplacian(kf2, T2)
      - fvm::laplacian(k2, T2)
      == -Q12
    );
    EEqn2.solve();

    //Generate mixture temperature
    T = (rhocpalpha1*T1 + rhocpalpha2*T2)/(rhocpalpha1 + rhocpalpha2); 

/*   
    const volScalarField DeltaT1 = mag(Q12)/rhocpalpha1;
    const volScalarField DeltaT2 = mag(Q12)/rhocpalpha2;
    const volScalarField DeltaTRatio = max( DeltaT1/mag(T2-T1), DeltaT2/mag(T1-T2) );
    const scalar RelaxationFactor = 1.0 / gMax( DeltaTRatio );
    Info<< endl << "Required relaxation factor: " << RelaxationFactor << endl;

    Info<< "Min/Max T1: " << gMin(T1) << "    " << gMax(T1) << endl;
    Info<< "Min/Max T2: " << gMin(T2) << "    " << gMax(T2) << endl;
    Info<< "Min/Max T: " << gMin(T) << "    " << gMax(T) << endl;
*/    
       
    
