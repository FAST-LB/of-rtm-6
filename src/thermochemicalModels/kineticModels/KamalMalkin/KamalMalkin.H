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
    Foam::KamalMalkin

Description
    Alexander Bernath, KIT-FAST-LBT, June 2017
    alexander.bernath@kit.edu

SourceFiles
    KamalMalkin.C
    KamalMalkinTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef KamalMalkin_H
#define KamalMalkin_H

#include "kineticModel.H"
#include "dimensionedTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticModels
{

/*---------------------------------------------------------------------------*\
                      Class KamalMalkin Declaration
\*---------------------------------------------------------------------------*/

class KamalMalkin
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

        //- Model coeffient n (user-supplied) [-]
        dimensionedScalar n_;

        //- Model coeffient m (user-supplied) [-]
        dimensionedScalar m_;

//        //- Name of density field
//        word rhoName_;
//
//        //- Name of dynamic viscosity field
//        word muName_;
//
//        //- Name of kinematic viscosity field
//        word nuName_;


    // Private Member Functions

//        //- Apply
//        template<class RhoFieldType>
//        void apply
//        (
//            scalarField& Udiag,
//            vectorField& Usource,
//            const scalarField& V,
//            const RhoFieldType& rho,
//            const scalarField& mu,
//            const vectorField& U
//        ) const;
//
//        //- Apply
//        template<class RhoFieldType>
//        void apply
//        (
//            tensorField& AU,
//            const RhoFieldType& rho,
//            const scalarField& mu,
//            const vectorField& U
//        ) const;

        //- Disallow default bitwise copy construct
        KamalMalkin(const KamalMalkin&);

        //- Disallow default bitwise assignment
        void operator=(const KamalMalkin&);


public:

    //- Runtime type information
    TypeName("KamalMalkin");

    //- Constructor
    KamalMalkin
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~KamalMalkin();


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
//    #include "KamalMalkinTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
