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
    Foam::regressionModel

Description
    Alexander Bernath, KIT-FAST-LBT, June 2017
    alexander.bernath@kit.edu

SourceFiles
    regressionModel.C

\*---------------------------------------------------------------------------*/

#ifndef regressionModel_H
#define regressionModel_H

#include "TgModel.H"
#include "dimensionedTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TgModels
{

/*---------------------------------------------------------------------------*\
                      Class regressionModel Declaration
\*---------------------------------------------------------------------------*/

class regressionModel
:
    public TgModel
{
    // Private data

    	//- Model coeffient Tg0 -> Tg at zero cure (user-supplied) [K]
        dimensionedScalar Tg0_;

        //- Model coeffient G1 (user-supplied) [-]
        dimensionedScalar G1_;

        //- Model coeffient G1 (user-supplied) [-]
        dimensionedScalar G2_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regressionModel(const regressionModel&);

        //- Disallow default bitwise assignment
        void operator=(const regressionModel&);


public:

    //- Runtime type information
    TypeName("regressionModel");

    //- Constructor
    regressionModel
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~regressionModel();


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

//#ifdef NoRepository
//    #include "regressionModelTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
