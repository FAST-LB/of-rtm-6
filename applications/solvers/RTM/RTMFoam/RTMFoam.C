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
    RTMFoam

    compressibleInterFoam extended by reaction kinetics, viscosity modeling
    and drag due to fibers.

    Alexander Bernath - KIT-FAST-LBT
    alexander.bernath@kit.edu
    
    Julian Seuffert - KIT-FAST-LBT
    julian.seuffert@kit.edu

Description
    CRTM solver for 2 compressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.
    
    Anisotropic porous media for multiple zones is included,
    Kinetic modeling for local viscosity is added,
    NOTE volume averaged velocity is treated by adjusting phi directly.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "compressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "IOkineticModelList.H"
#include "IOTgModelList.H"
#include "IOporosityModelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createZones.H" 
    
    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();
    const volScalarField& psi1 = mixture.thermo1().psi();
    const volScalarField& psi2 = mixture.thermo2().psi();    
    
    #include "createKineticFields.H"
    
    if (!LTS)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        { 	              
           
            // Correct phi with fiberVolFraction
            phi /= (scalar(1) - fvc::interpolate(fiberVolFraction));	

            #include "alphaControls.H"
            #include "porousCompressibleAlphaEqnSubCycle.H"

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
            
            mu = alpha1*mixture.thermo1().mu() + alpha2*mixture.thermo2().mu();

            // --- Pressure corrector loop
            while (pimple.correct())
            {                
                #include "pEqn.H"
                p.correctBoundaryConditions();
            }            
            
            Info<< "max(U) " << max(mag(U)).value() << endl;
            Info<< "min(p_rgh): " << min(p_rgh).value() << "  max(p_rgh): " << max(p_rgh).value() << endl;

            if (pimple.turbCorr())
            {
                turbulence.correct();
            }
        }
            
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
        
        #include "calcFillingTime.H"    
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
