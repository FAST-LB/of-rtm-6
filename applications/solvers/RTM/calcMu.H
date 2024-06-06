
if (calcMu)
{    
    fvScalarMatrix resinAgeEqn
    (
        fvm::ddt(rho, resinAge)
      + fvm::div(rhoPhi, resinAge)
      - fvm::Sp(contErr, resinAge)
    );

    resinAgeEqn.relax();
    resinAgeEqn.solve();
    
    resinAge += runTime.deltaT();
    
    Polynomial<8> muCoeffs(mixture.thermo1().lookupType<Polynomial<8>>("muCoeffs"));
    
    forAll(mesh.C(), cellI)  
    {
        dimensionedScalar mu1("mu1", dimMass/dimLength/dimTime, muCoeffs.value(resinAge[cellI]));
        mu[cellI] = alpha1[cellI]*mu1.value();
    }
    
    mu.ref() += alpha2*mixture.thermo2().mu();
    Info<< "min/max mu: " << min(U.db().lookupObject<volScalarField>("muFluid").internalField()).value()
                    << "/" << max(U.db().lookupObject<volScalarField>("muFluid").internalField()).value() << endl;         
}
else
{
    mu = alpha1*mixture.thermo1().mu() + alpha2*mixture.thermo2().mu();
}
