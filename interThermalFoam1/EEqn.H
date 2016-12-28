    //Get conductivity on cell faces
    const surfaceScalarField alpha1f = fvc::interpolate(alpha1);
    const surfaceScalarField kf = alpha1f*k1 + (1.0-alpha1f)*k2;

    fvScalarMatrix EEqn
    (
        fvm::ddt(rhocp, T)
      + fvm::div(rhocpPhi, T)
     ==
        fvm::laplacian(kf, T)
    );

    EEqn.solve();
