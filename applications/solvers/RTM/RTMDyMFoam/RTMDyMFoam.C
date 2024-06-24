/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    RTMDyMFoam

    compressibleInterDyMFoam extended by reaction kinetics, viscosity modeling
    and drag due to fibers.

    Alexander Bernath - KIT-FAST-LBT
    alexander.bernath@kit.edu
    
    Julian Seuffert - KIT-FAST-LBT
    julian.seuffert@kit.edu

Description
    CRTM solver for 2 compressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.
    
    Anisotropic porous media for multiple zones is included,
    Kinetic modeling for local viscosity is added,
    Fiber volume fraction is updated each time step depending on cell volume change,
    NOTE volume averaged velocity is treated by adjusting phi directly

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "compressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "IOkineticModelList.H"
#include "IOTgModelList.H"
#include "IOporosityModelList.H"
#include "displacementMotionSolver.H"
#include "volPointInterpolation.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"    
    #include "createFields.H"
    #include "createZones.H"

    //#include "createForces.H"
    
    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();
    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();    
    
    #include "createKineticFields.H"
    #include "createUf.H"
    #include "CourantNo.H"  
    #include "initContinuityErrs.H"
    #include "setInitialDeltaT.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {        
        #include "readDyMControls.H"

        // Store divU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        volScalarField divU("divU0", fvc::div(fvc::absolute(phi, U)));

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "deformationCourantNo.H"
            #include "setDeformationDeltaT.H"
        }
        
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

	// --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        { 
            /*
            if (pimple.firstIter())
            {
                #include "p_rghPredictor.H"
            }
            */
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                //#include "p_rghPredictor.H"
                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                mesh.update();

                if (mesh.changing())
                {
                        
                    MRF.update();

                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;
                            
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf;

                        #include "correctPhi_rgh.H"

                        // Make the fluxes relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    #include "updatePorousZones.H"
                }
            }
            
            // Correct phi with fiberVolFraction
            phi /= (scalar(1) - fvc::interpolate(fiberVolFraction));	

            #include "alphaControls.H"
            #include "compressibleAlphaEqnSubCycle.H"

            turbulence.correctPhasePhi();           

            #include "porousUEqn.H"

            {
                mixture.correctThermo();
                mixture.correct();
            }
            
            if (kineticZones.active())
            {  
                #include "cureEqn.H"
            }
            
            #include "calcMu.H"
            
            // --- Pressure corrector loop
            while (pimple.correct() )
            {
                if (mesh.changing())
                {
                    #include "pEqnDyM.H"
                }
                else
                {
                    #include "pEqn.H"
                }
                p.correctBoundaryConditions();
            }  
            
            
            Info<< "max(U) " << max(mag(U)).value() << endl;
            Info<< "min(p_rgh): " << min(p_rgh).value() << "  max(p_rgh): " << max(p_rgh).value() << endl;

            if (pimple.turbCorr())
            {
                turbulence.correct();
            }            
        }       
        
        #include "calcFillingTime.H"

        if (kineticZones.active())
        {
            kineticZones.calcCureRate(cure,T,Tg,cureRate);                    
            kineticZones.calcCure(cure,cureRate);
                    
            forAll(alpha1, cellI)
            {
                if (alpha1[cellI] > 0.001)
                {
                    cureEff[cellI] = cure[cellI];
                } 
                else 
                {
                    cureEff[cellI] = 0.0;
                }
            }

            if (TgZones.active())
            {
                TgZones.calcTg(cure,T,Tg);
                forAll(alpha1, cellI)
                {
                    if (alpha1[cellI] > 0.001)
                    {
                        TgEff[cellI] = Tg[cellI];
                    } 
                    else 
                    {
                        TgEff[cellI] = 0.0;
                    }
                }
                kineticZones.calcMaterialState(cureEff,T,TgEff,materialState);
            } 
            else
            {
                kineticZones.calcMaterialState(cureEff,materialState);
            } 
            Info<< "min/max cureEff: " << min(cureEff.internalField()).value()
                    << "/" << max(cureEff.internalField()).value() << endl;
            Info<< "min/max TgEff: " << min(TgEff.internalField()).value()
                    << "/" << max(TgEff.internalField()).value() << endl;
            Info<< "min/max mu: " << min(U.db().lookupObject<volScalarField>("thermo:mu").internalField()).value()
                    << "/" << max(U.db().lookupObject<volScalarField>("thermo:mu").internalField()).value() << endl;            
            Info<< "min/max mu.resin: " << min(U.db().lookupObject<volScalarField>("thermo:mu.resin").internalField()).value()
                    << "/" << max(U.db().lookupObject<volScalarField>("thermo:mu.resin").internalField()).value() << endl;           
        }
        
        runTime.write();
                
        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;            
        
        #include "endSimulation.H"            
    }
        
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
