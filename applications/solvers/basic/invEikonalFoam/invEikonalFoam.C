/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    invEikonalFoam

Description
    Solves the inverse wall distance equation, e.g. for fast (mold)filling approximations.
    To obtain the distance field, postprocessing is mandatory.

Reference
    Fares and Schroeder. A differential equation for approximate wall distance.
    International Journal for Numerical Methods in Fluids. 2002. 39(8):743-764.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "fvOptions.H"
#include "simpleControl.H"
//#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating inverse wall distance\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            U = fvc::grad(G);
            
            fvScalarMatrix GEqn
            (
                tau * fvc::ddt(G)
                + fvm::laplacian(G, G)
                + (sigma-1.0)*G*0.5*(fvm::laplacian(G)+fvc::div(U))
             ==
                  (1.0 + 2.0*sigma)*pow(h/H,2.0)*pow(G,4.0)
            );
            
            GEqn.solve();
            G.relax();
        }

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
