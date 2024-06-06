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
    Foam::Grindling

Description
    Alexander Bernath, KIT-FAST-LBT, June 2017
    alexander.bernath@kit.edu

SourceFiles
    Grindling.C

\*---------------------------------------------------------------------------*/

#ifndef Grindling_H
#define Grindling_H

#include "kineticModel.H"
#include "dimensionedTensor.H"
#include "DimensionedScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticModels
{

/*---------------------------------------------------------------------------*\
                      Class Grindling Declaration
\*---------------------------------------------------------------------------*/

class Grindling
:
    public kineticModel
{
    // Private data

        //- universal gas constant [J/(mol*K)]
        dimensionedScalar R_;

        //- Model coeffient A1 (user-supplied) [1/s]
        dimensionedScalar A1_;

        //- Model coeffient A2 (user-supplied) [1/s]
        dimensionedScalar A2_;

        //- Model coeffient E1 (user-supplied) [J/mol]
        dimensionedScalar E1_;

        //- Model coeffient E2 (user-supplied) [J/mol]
        dimensionedScalar E2_;

        //- Model coeffient n1 (user-supplied) [-]
        dimensionedScalar n1_;

        //- Model coeffient n2 (user-supplied) [-]
        dimensionedScalar n2_;

        //- Model coeffient m (user-supplied) [-]
        dimensionedScalar m_;

        //- Model coeffient m (user-supplied) [1/s]
        dimensionedScalar K2dTg_;

        //- Model coeffient c1 (user-supplied) [-]
        dimensionedScalar c1_;

        //- Model coeffient c2 (user-supplied) [K]
        dimensionedScalar c2_;

        //- Model coeffient deltaTg (user-supplied) [K]
        dimensionedScalar deltaTg_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Grindling(const Grindling&);

        //- Disallow default bitwise assignment
        void operator=(const Grindling&);


public:

    //- Runtime type information
    TypeName("Grindling");

    //- Constructor
    Grindling
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~Grindling();


    // Member Functions

        //- Calculate the cure rate
        virtual void calcCureRate
        (
            const volScalarField& cure,
            const volScalarField& T,
            const volScalarField& Tg,
			volScalarField& cureRate
        ) const;

    // I-O

        //- Write
        bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "GrindlingTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
