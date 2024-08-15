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

#include "Grindling.H"

#include "addToRunTimeSelectionTable.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace kineticModels
    {
        defineTypeNameAndDebug(Grindling, 0);
        addToRunTimeSelectionTable(kineticModel, Grindling, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticModels::Grindling::Grindling
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
    n1_("n1", dimless, coeffs_),
    n2_("n2", dimless, coeffs_),
    m_("m", dimless, coeffs_),
    K2dTg_("K2dTg", dimless/dimTime, coeffs_),
    c1_("c1", dimless, coeffs_),
    c2_("c2", dimTemperature, coeffs_),
    deltaTg_("deltaTg", dimTemperature, coeffs_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticModels::Grindling::~Grindling()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticModels::Grindling::calcCureRate
(
    const volScalarField& cure,
    const volScalarField& T,
    const volScalarField& Tg,
	volScalarField& cureRate
) const
{
	volScalarField K1 = A1_*exp(-E1_/(R_*T));
	volScalarField K2 = A2_*exp(-E2_/(R_*T));

	volScalarField E1diff = R_*pow(Tg,scalar(2.0))*c1_/c2_;
	volScalarField E2diff = R_*c1_*c2_*pow(Tg+deltaTg_,scalar(2.0))/pow(c2_+deltaTg_,scalar(2.0));

	volScalarField K2diff = K2dTg_*exp(c1_*(T-Tg)/(c2_+T-Tg));

	volScalarField Keff = scalar(1.0)/(scalar(1.0)/K2diff+scalar(1.0)/K2);

	cureRate = K1*pow(max(scalar(1.0)-cure,scalar(0.0)),n1_)+Keff*pow(max(scalar(0.0),cure),m_)*pow(max(scalar(1.0)-cure,scalar(0)),n2_);

//    return cureRate;

//	forAll(K2diff,cellI)
//	{
//        if (T[cellI]>(Tg[cellI]+deltaTg_.value())) {
//            K2diff[cellI] = ((K2dTg_*exp(c1_*deltaTg_/(c2_+deltaTg_)))*exp(-E2diff[cellI]/R_.value()*(scalar(1.0)/T[cellI]-scalar(1.0)/(Tg[cellI]+deltaTg_.value())))).value();
//        } else if (T[cellI]<Tg[cellI]) {
//            K2diff[cellI] = (K2dTg_*exp(-E1diff[cellI]/R_.value()*(scalar(1.0)/T[cellI]-scalar(1.0)/Tg[cellI]))).value();
//        } else {
//            K2diff[cellI] = (K2dTg_*exp(c1_*(T[cellI]-Tg[cellI])/(c2_+T[cellI]-Tg[cellI]))).value();
//        }
//        Keff[cellI] = scalar(1.0)/(scalar(1.0)/K2diff[cellI]+scalar(1.0)/K2[cellI]);
//        cureRate[cellI] = (K1[cellI]*pow(max(scalar(1.0)-cure[cellI],scalar(0.0)),n1_)
//        		        + Keff[cellI]*pow(cure[cellI],m_)*pow(max(scalar(1.0)-cure[cellI],scalar(0)),n2_)).value();
//	}

//	scalarField Keff = scalar(1.0)/(scalar(1.0)/K2diff+scalar(1.0)/K2);
//
//	cureRate = K1*pow(max(scalar(1.0)-cure,scalar(0.0)),n1_)+Keff*pow(cure,m_)*pow(max(scalar(1.0)-cure,scalar(0)),n2_);
}

bool Foam::kineticModels::Grindling::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
