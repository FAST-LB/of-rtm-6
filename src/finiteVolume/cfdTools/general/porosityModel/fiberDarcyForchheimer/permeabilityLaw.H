/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

Class
    Foam::permeabilityLaw

Description
    permeability model base class


\*---------------------------------------------------------------------------*/

#ifndef permeabilityLaw_H
#define permeabilityLaw_H

#include "fvMesh.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "porosityModel.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

namespace Foam
{
    
class permeabilityLaw  
{       
   // Private data
        //- Name
        const fvMesh& mesh_;
        
        //- Name
        const word name_;        
        
        //- cellZoneID
        const label& cellZoneID_;
        
        //- Disallow default bitwise copy construct
        permeabilityLaw(const permeabilityLaw&);

        //- Disallow default bitwise assignment
        void operator=(const permeabilityLaw&);
        
protected:

        const dictionary& dict_;
        
    // Private Member Functions

public:

    //- Runtime type information
    TypeName("permeabilityLaw");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            permeabilityLaw,
            dictionary,
            (
                const fvMesh& mesh,
                const word name,
                const label& cellZoneID,
                const dictionary& dict
            ),
            (mesh, name, cellZoneID, dict)
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<permeabilityLaw> New
        (
            const fvMesh& mesh,
            const word& name,
            const label& cellZoneID,
            const dictionary& dict
        );

    
    //- Constructor
        
        permeabilityLaw
        ( 
            const fvMesh& mesh,
            const word& name,
            const label& cellZoneID,
            const dictionary& dict
        );

    //- Destructor
        
        virtual ~permeabilityLaw()
        {}
    
    // Member Functions
    
        //- Return mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        } 

        //- Return name
        const word& name() const
        {
            return name_;
        } 
        
        //- Return name
        const label& cellZoneID() const
        {
            return cellZoneID_;
        } 
        
                
        //- Darcy
        virtual tensor D(const scalar& FVF, const label& j) = 0;    
        
        //- Reference to Darcy component
        virtual tmp<tensorField> D(const scalarField& FVF) = 0; 
        
        virtual tmp<tensorField> D() = 0;
        
        //- Forchheimer
        virtual tensor F(const scalar& FVF, const label& j) = 0; 
        
        //- Reference to Forchheimer component
        virtual tmp<tensorField> F(const scalarField& FVF) = 0;
        
        virtual tmp<tensorField> F() = 0;   
    
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //    
