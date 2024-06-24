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

Global
    deformationCourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/
bool deformationCo = false;
runTime.controlDict().readIfPresent("deformationCo", deformationCo);

scalar maxDeformationCo
(
    runTime.controlDict().lookupOrDefault("maxDeformationCo", 1.0)
); 
    
scalar deformationCoNum = 0.0;
scalar meanDeformationCoNum = 0.0;

if (mesh.nInternalFaces() && runTime.timeIndex() > 1 && deformationCo == true)
{
    volScalarField& rhoSolid
    ( 
        const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("rhoSolid"))
    );   
    
    // Get needed reference to the mesh
    const edgeList& edges = mesh.edges();
    const pointField& pp = mesh.points();
    const labelListList& cEdges = mesh.cellEdges();

    // Make a list with all the edge lenghts
    scalarField eLengths(edges.size(), 0.0);

    forAll (edges, edgei)
    {
        eLengths[edgei] = edges[edgei].mag(pp);
    }

    // Find the minimum per cell
    scalarField minLengthCells(cEdges.size(), GREAT);
    scalarField Cd(cEdges.size(), 0.0);
    scalarField dT(cEdges.size(), 0.0);
    scalarField deformationCoNumField(cEdges.size(), 0.0);
    
    forAll (cEdges, celli)
    {
        const labelList& ce = cEdges[celli];

        forAll (ce, edgei)
        {
            minLengthCells[celli] = Foam::min(minLengthCells[celli], eLengths[ce[edgei]]);
        }
        
        dimensionedScalar twoMuLambdaCell = mesh.lookupObject<volScalarField>("twoMuLambda")[celli];
        
        dimensionedScalar rhoSolidCell = rhoSolid[celli];
        
        Cd[celli] = sqrt(twoMuLambdaCell/rhoSolidCell).value();
        
        dT[celli] = minLengthCells[celli]/Cd[celli];
        
        deformationCoNumField[celli] = runTime.deltaTValue()/dT[celli];
    }   
    
    //Local mass scaling
    if (deformationMassScaling)
    {    
        scalar maxMassScaling
        (
            readScalar(runTime.controlDict().lookup("maxMassScaling"))
        );
        
        volScalarField scaleFactor0(scaleFactor);

        forAll(mesh.C(), celli)
        {            
            //scale with CoNum, maxMassscaling
            scaleFactor[celli] *= Foam::max(1.0/1.2, Foam::pow(deformationCoNumField[celli]/maxDeformationCo, 2))*1.2;
            scaleFactor[celli] = Foam::min(maxMassScaling, scaleFactor[celli]);
            rhoSolid[celli] *= scaleFactor[celli]/scaleFactor0[celli];
                
            //adjust CoNum field
            deformationCoNumField[celli] /= Foam::sqrt(scaleFactor[celli]/scaleFactor0[celli]);                           
        }                    
                
        Info << " Maximum mass scaling: " << max(scaleFactor).value() << endl;        
    }
    
    meanDeformationCoNum =
                average(deformationCoNumField);
    deformationCoNum = max(deformationCoNumField);
    
    Info<< "Deformation Courant Number mean: " << meanDeformationCoNum
    << " max: " << deformationCoNum << endl;
}



// ************************************************************************* //
