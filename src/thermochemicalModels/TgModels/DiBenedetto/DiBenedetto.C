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

#include "addToRunTimeSelectionTable.H"
#include "DiBenedetto.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace TgModels
    {
        defineTypeNameAndDebug(DiBenedetto, 0);
        addToRunTimeSelectionTable(TgModel, DiBenedetto, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TgModels::DiBenedetto::DiBenedetto
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    TgModel(name, modelType, mesh, dict, cellZoneName),
    Tg0_("Tg0", dimTemperature, coeffs_),
    TgInf_("TgInf", dimTemperature, coeffs_),
	lambda_("lambda", dimless, coeffs_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TgModels::DiBenedetto::~DiBenedetto()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TgModels::DiBenedetto::calcTg
(
    const volScalarField& cure,
    const volScalarField& T,
    volScalarField& Tg
) const
{
	Tg = Tg0_+(TgInf_-Tg0_)*lambda_*cure/(scalar(1)-(scalar(1)-lambda_)*cure);
}

bool Foam::TgModels::DiBenedetto::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
