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
    HeleShaw

Description
    HeleShaw effective permeability with K = h^2 / 12

SourceFiles
    HeleShaw.C

\*---------------------------------------------------------------------------*/

#ifndef HeleShaw_H
#define HeleShaw_H

#include "permeabilityLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class HeleShaw Declaration
\*---------------------------------------------------------------------------*/

class HeleShaw
:
    public permeabilityLaw
{
    // Private data
    
        //- Name
        const fvMesh& mesh_;
        
        //- Dictionary        
        const dictionary& HeleShawDict;
        
        //- cellZoneID
        const label cellZoneID_;
       
        scalar H0_;
        scalar Hdot_;        
        
        //- viscosity (mu in thermophysicalProperties should be set to 1 Pas)
        scalar mu_;
        bool localGapHeight_;
        
        List<word> patchNames_;
        List<word> faceZoneNames_;
        List<label> patchIDs_;
        List<label> faces1_;
        List<label> faces2_;
        
        scalarField height_;
        scalar timeValue_;
        
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        HeleShaw(const HeleShaw&);

        //- Disallow default bitwise assignment
        void operator=(const HeleShaw&);

public:

    //- Runtime type information
    TypeName("HeleShaw");

    // Static data members


    // Constructors

        //- Construct from dictionary
        HeleShaw
        (
            const fvMesh& mesh,
            const word& name,
            const label& cellZoneID,
            const dictionary& dict
        );

    // Destructor

        virtual ~HeleShaw();


    // Member Functions
        
        //- Return height scalar
        virtual scalar H(const label& j);
        
        //- Return height scalar
        virtual scalar H(const scalar& t, const label& j);
        
        //- Return height scalarField
        virtual scalarField H(const scalar& t);
        
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
