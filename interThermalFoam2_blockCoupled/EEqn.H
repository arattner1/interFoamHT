    //Get conductivity on cell faces
    const surfaceScalarField alpha1f = fvc::interpolate(alpha1);
    const surfaceScalarField kf1 = alpha1f*k1;
    const surfaceScalarField kf2 = (1.0-alpha1f)*k2;

    //Limiting to keep problem conditioned, otherwise the equations get a 
    //divided by 0 type issue
    const scalar eps = 1E-8;
    rhocpalpha1 = max( rhocpalpha1, eps*rho1*cp1);
    rhocpalpha2 = max( rhocpalpha2, eps*rho2*cp2);


    //Interphase transport model
    //Interface area density (m^2/m^3)
    const volVectorField gradAlpha = fvc::grad(alpha1);
    const volScalarField Aint = mag(gradAlpha);
    //Get normal to interface
    const dimensionedScalar deltaN("deltaN", dimensionSet(0,-1,0,0,0,0,0), 1E-16);
	const volVectorField ns = gradAlpha/(mag(gradAlpha) + deltaN);
	//Average k in each cell
	const volScalarField kc = (alpha1*k1 + (1.0-alpha1)*k2);
    //hA_effective in each cell
    const volScalarField hA_eff = 0.0 * kc*mag(gradAlpha)*Aint;

    //Try block-matrix solve:
    //- Prepare block system
    BlockLduMatrix<vector2> blockM(mesh);
    //- Grab block diagonal and set it to zero
    Field<tensor2>& d = blockM.diag().asSquare();
    d = tensor2::zero;
    //- Grab linear off-diagonal and set it to zero
    Field<vector2>& u = blockM.upper().asLinear();
    Field<vector2>& l = blockM.lower().asLinear();
    u = vector2::zero;
    l = vector2::zero;
    //- Create the source term
    vector2Field blockB(mesh.nCells(), vector2::zero);
    //- Create the working solution variable
    vector2Field blockX(mesh.nCells(), vector2::zero);


    fvScalarMatrix EEqn1
    (
        fvm::ddt(rhocpalpha1, T1)
      + fvm::div(rhocpPhi1, T1)
      - fvm::laplacian(kf1, T1)
      == hA_eff*T2 - fvm::Sp(hA_eff, T1)
    );

    fvScalarMatrix EEqn2
    (
        fvm::ddt(rhocpalpha2, T2)
      + fvm::div(rhocpPhi2, T2)
      - fvm::laplacian(kf2, T2)
      == hA_eff*T1 - fvm::Sp(hA_eff, T2)
    );

    //- Insert equations into block Matrix
    blockMatrixTools::insertEquation(0, EEqn1, blockM, blockX, blockB);
    blockMatrixTools::insertEquation(1, EEqn2, blockM, blockX, blockB);

    //- Add off-diagonal terms and remove from Block source
    forAll(d, i)
    {
        d[i](0,1) = -hA_eff.internalField()[i]*mesh.V()[i];
        d[i](1,0) = -hA_eff.internalField()[i]*mesh.V()[i];
        blockB[i][0] -= hA_eff.internalField()[i]*blockX[i][1]*mesh.V()[i];
        blockB[i][1] -= hA_eff.internalField()[i]*blockX[i][0]*mesh.V()[i];
    }

    //- Block coupled solver call
    BlockSolverPerformance<vector2> solverPerf =
        BlockLduSolver<vector2>::New
        (
            word("blockVar"),
            blockM,
            mesh.solutionDict().solver("T")
        )->solve(blockX, blockB);
    solverPerf.print();
        
    //- Retrieve solution
    blockMatrixTools::blockRetrieve(0, T1.internalField(), blockX);
    blockMatrixTools::blockRetrieve(1, T2.internalField(), blockX);
    T1.correctBoundaryConditions();
    T2.correctBoundaryConditions();


//    solve( EEqn1 == hA_eff*(T2 - T1) );
//    solve( EEqn2 == hA_eff*(T1 - T2) );
    
    
    
    
     
//    const dimensionedScalar deltat("dt", dimensionSet(0,0,1,0,0,0,0), runTime.deltaTValue() );
//    const volScalarField DeltaT1 = deltat*hA_eff*(T2-T1)/(rhocpalpha1);
//    const volScalarField DeltaT2 = deltat*hA_eff*(T1-T2)/(rhocpalpha2);    
//    const volScalarField DeltaTRatio = max( DeltaT1/(T2-T1), DeltaT2/(T1-T2) );
//    const scalar RelaxationFactor = 1.0 / gMax( DeltaTRatio );
//    Info<< endl << "Required relaxation factor: " << RelaxationFactor << endl << endl;
    Info<< endl << "Min/Max T1: " << gMin(T1) << "    " << gMax(T1) << endl;
    Info<< "Min/Max T2: " << gMin(T2) << "    " << gMax(T2) << endl;

    
    T = (rhocpalpha1*T1 + rhocpalpha2*T2)/(rhocpalpha1 + rhocpalpha2);            
    
//    T1 = pos(alpha1 - eps)*T1 + neg(alpha1 - eps)*T2;
//    T2 = pos(alpha2 - eps)*T2 + neg(alpha2 - eps)*T1;
