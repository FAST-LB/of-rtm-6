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
    Gebart

Description
    Linear elastic permeability

SourceFiles
    Gebart.C

\*---------------------------------------------------------------------------*/

#ifndef Gebart_H
#define Gebart_H

#include "permeabilityLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Gebart Declaration
\*---------------------------------------------------------------------------*/

class Gebart
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

        //-fiber volume fraction for the specified permeability data
        //Only used for local treatment of fiber volume fraction
        scalar vf0_;

        //-	maximum value for the fiber volume fraction
        //	needed for the calculation of the local fiber volume fraction as proposed by Gebart (1992)
        scalar vfmax_;
        
        scalar constDA_;
        scalar constDB1_;
        scalar constDB2_;
        
        scalar constFA_;
        scalar constFB1_;
        scalar constFB2_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        Gebart(const Gebart&);

        //- Disallow default bitwise assignment
        void operator=(const Gebart&);

public:

    //- Runtime type information
    TypeName("Gebart");

    // Static data members


    // Constructors

        //- Construct from dictionary
        Gebart
        (
            const fvMesh& mesh,
            const word& name,
            const label& cellZoneID,
            const dictionary& dict
        );


    // Destructor

        virtual ~Gebart();


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
