    
    fvVectorMatrix UEqn
    (
        fvm::div(fvc::interpolate(rho)*(mesh.Sf()&Uf), U)      
    );  

    // Include the porous media resistance and solve the momentum equation
    // either implicit in the tensorial resistance or transport using by
    // including the spherical part of the resistance in the momentum diagonal
    pZones.addResistance(UEqn);
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    
    solve(
        - fvm::laplacian(rAUf, p_rgh)        
    );
    
    p = p_rgh + (alpha1*rho1 + alpha2*rho2)*gh;
    p_rgh = p - (alpha1*rho1 + alpha2*rho2)*gh;
    
    Info<< "min(p_rgh): " << min(p_rgh).value() << "  max(p_rgh): " << max(p_rgh).value() << endl; 

