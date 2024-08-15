/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Foam::DiBenedetto

Description
    Alexander Bernath, KIT-FAST-LBT, June 2017
    alexander.bernath@kit.edu

SourceFiles
    DiBenedetto.C

\*---------------------------------------------------------------------------*/

#ifndef DiBenedetto_H
#define DiBenedetto_H

#include "TgModel.H"
#include "dimensionedTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TgModels
{

/*---------------------------------------------------------------------------*\
                      Class DiBenedetto Declaration
\*---------------------------------------------------------------------------*/

class DiBenedetto
:
    public TgModel
{
    // Private data

    	//- Model coeffient Tg0 -> Tg at zero cure (user-supplied) [K]
        dimensionedScalar Tg0_;

        //- Model coeffient TgInf -> Tg at full cure (user-supplied) [K]
        dimensionedScalar TgInf_;

        //- Model coeffient lambda (user-supplied) [-]
        dimensionedScalar lambda_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        DiBenedetto(const DiBenedetto&);

        //- Disallow default bitwise assignment
        void operator=(const DiBenedetto&);


public:

    //- Runtime type information
    TypeName("DiBenedetto");

    //- Constructor
    DiBenedetto
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~DiBenedetto();


    // Member Functions

        //- Calculate Tg from model
        virtual void calcTg
        (
            const volScalarField& cure,
            const volScalarField& T,
            volScalarField& Tg
        ) const;

    // I-O

        //- Write
        bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TgModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
