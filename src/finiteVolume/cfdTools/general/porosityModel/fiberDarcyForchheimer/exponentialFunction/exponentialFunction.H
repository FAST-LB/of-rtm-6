/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    exponentialFunction

Description
    Linear elastic permeability

SourceFiles
    exponentialFunction.C

\*---------------------------------------------------------------------------*/

#ifndef exponentialFunction_H
#define exponentialFunction_H

#include "permeabilityLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class exponentialFunction Declaration
\*---------------------------------------------------------------------------*/

class exponentialFunction
:
    public permeabilityLaw
{
    // Private data
    
        //- Name
        const fvMesh& mesh_;
        
        //- cellZoneID
        const label cellZoneID_;

        
        //- Darcy coefficient
        dimensionedTensor DA_;

        //- Forchheimer coefficient
        dimensionedTensor FO_;       

         //- exponentialFunction coefficient C0
        scalar C01_;
        scalar C02_;
        scalar C03_;
        //- exponentialFunction coefficient C1
        scalar C11_;
        scalar C12_;
        scalar C13_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        exponentialFunction(const exponentialFunction&);

        //- Disallow default bitwise assignment
        void operator=(const exponentialFunction&);

public:

    //- Runtime type information
    TypeName("exponentialFunction");

    // Static data members


    // Constructors

        //- Construct from dictionary
        exponentialFunction
        (
            const fvMesh& mesh,
            const word& name,
            const label& cellZoneID,
            const dictionary& dict
        );


    // Destructor

        virtual ~exponentialFunction();


    // Member Functions
        
        //- Return yield stress given cellID
        virtual tensor D(const scalar& FVF, const label& j);    

        //- Reference to Darcy component
        virtual tmp<tensorField> D(const scalarField& FVF); 
        
        virtual tmp<tensorField> D();
        
        //- Return yield stress given cellID
        virtual tensor F(const scalar& FVF, const label& j);    
        
        //- Reference to Forchheimer component
        virtual tmp<tensorField> F(const scalarField& FVF);
        
        virtual tmp<tensorField> F();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
