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
    tabularData

Description
    Linear elastic permeability

SourceFiles
    tabularData.C

\*---------------------------------------------------------------------------*/

#ifndef tabularData_H
#define tabularData_H

#include "permeabilityLaw.H"
#include "spline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tabularData Declaration
\*---------------------------------------------------------------------------*/

class tabularData
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

        // Values needed for the calculation of the interpolation of measured data
	scalarField messfvg_;
	scalarField messk1_;
	scalarField messk2_;
	scalarField messk3_;
	scalarField A1_, A2_, A3_, B1_, B2_, B3_;
        
        //Values for spline interpolation
	tk::spline KK1;
	tk::spline KK2;
	tk::spline KK3;

	word interpolationType_;            
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        tabularData(const tabularData&);

        //- Disallow default bitwise assignment
        void operator=(const tabularData&);

public:

    //- Runtime type information
    TypeName("tabularData");

    // Static data members


    // Constructors

        //- Construct from dictionary
        tabularData
        (
            const fvMesh& mesh,
            const word& name,
            const label& cellZoneID,
            const dictionary& dict
        );


    // Destructor

        virtual ~tabularData();


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
