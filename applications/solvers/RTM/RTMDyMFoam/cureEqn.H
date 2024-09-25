{
    fvScalarMatrix cureEqn
    (
        fvm::ddt(rho, cure)
      + fvm::div(rhoPhi, cure)
      - fvm::Sp(contErr, cure)
//      - fvm::laplacian(mixture.alphaEff(turbulence->mut()), cure)
//      + (
//            fvc::div(fvc::absolute(phi, U), p)
//          + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
//        )
//       *(
//           alpha1/mixture.thermo1().Cv()
//         + alpha2/mixture.thermo2().Cv()
//        )
    );

    cureEqn.relax();
    cureEqn.solve();
    
    Info<< "min(cure)/max(cure): " << min(cure).value() << "/" << max(cure).value() << endl;
}
