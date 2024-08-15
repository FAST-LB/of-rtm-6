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
#include "KamalMalkin.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace kineticModels
    {
        defineTypeNameAndDebug(KamalMalkin, 0);
        addToRunTimeSelectionTable(kineticModel, KamalMalkin, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticModels::KamalMalkin::KamalMalkin
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    kineticModel(name, modelType, mesh, dict, cellZoneName),
    R_("R", dimEnergy/dimMoles/dimTemperature, coeffs_),
    A1_("A1", dimless/dimTime, coeffs_),
	A2_("A2", dimless/dimTime, coeffs_),
    E1_("E1", dimEnergy/dimMoles, coeffs_),
	E2_("E2", dimEnergy/dimMoles, coeffs_),
    n_("n", dimless, coeffs_),
    m_("m", dimless, coeffs_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticModels::KamalMalkin::~KamalMalkin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticModels::KamalMalkin::calcCureRate
(
    const volScalarField& cure,
    const volScalarField& T,
    const volScalarField& Tg,
	volScalarField& cureRate
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
        forAll(cells,j)
        {
            const dimensionedScalar TC ("T", dimTemperature, T[cells[j]]);
            const dimensionedScalar K1 = A1_*exp(-E1_/(R_*TC));
            const dimensionedScalar K2 = A2_*exp(-E2_/(R_*TC));

            // Set the internal field values.
            cureRate[cells[j]] = ((K1+K2*pow(max(cure[cells[j]], scalar(0.0)), m_))*pow(max((scalar(1.0)-cure[cells[j]]), scalar(0)), n_)).value();
        }
    }
}

bool Foam::kineticModels::KamalMalkin::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
