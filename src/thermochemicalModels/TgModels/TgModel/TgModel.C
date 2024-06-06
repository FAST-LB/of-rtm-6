/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "TgModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TgModel, 0);
    defineRunTimeSelectionTable(TgModel, mesh);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

//Foam::label Foam::TgModel::fieldIndex(const label i) const
//{
//    label index = 0;
//    if (!coordSys_.R().uniform())
//    {
//        index = i;
//    }
//    return index;
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TgModel::TgModel
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.subDict(modelType + "Coeffs")),
    active_(true),
    zoneName_(cellZoneName),
    cellZoneIDs_()
{
    if (zoneName_ == word::null)
    {
        dict.lookup("active") >> active_;
        dict_.lookup("cellZone") >> zoneName_;
    }

    cellZoneIDs_ = mesh_.cellZones().findIndices(zoneName_);

    Info<< "    creating Tg model: " << name_ << endl;

    bool foundZone = !cellZoneIDs_.empty();
    reduce(foundZone, orOp<bool>());

    if (!foundZone && Pstream::master())
    {
        FatalErrorInFunction
            << "cannot find Tg cellZone " << zoneName_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TgModel::~TgModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::TgModel::writeData(Ostream& os) const
{
    return true;
}


bool Foam::TgModel::read(const dictionary& dict)
{
    active_ = readBool(dict.lookup("active"));
    coeffs_ = dict.subDict(type() + "Coeffs");
    dict.lookup("cellZone") >> zoneName_;
    cellZoneIDs_ = mesh_.cellZones().findIndices(zoneName_);

    return true;
}


// ************************************************************************* //
