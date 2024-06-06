/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "kineticModelList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticModelList::kineticModelList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<kineticModel>(),
    mesh_(mesh)
{
    reset(dict);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticModelList::~kineticModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::kineticModelList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "No kinetic models active" << endl;
    }

    return a;
}


void Foam::kineticModelList::reset(const dictionary& dict)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& modelDict = iter().dict();

            this->set
            (
                i++,
                kineticModel::New(name, mesh_, modelDict)
            );
        }
    }
}


bool Foam::kineticModelList::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        kineticModel& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::kineticModelList::writeData(Ostream& os) const
{
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    return os.good();
}

void Foam::kineticModelList::calcCure
(
    volScalarField& cure,
	const volScalarField& cureRate
)
{
    forAll(*this, i)
    {
        this->operator[](i).calcCure(cure,cureRate);
    }
}

void Foam::kineticModelList::calcCureRate
(
    const volScalarField& cure,
	const volScalarField& T,
	const volScalarField& Tg,
	volScalarField& cureRate
)
{
    forAll(*this, i)
    {
        this->operator[](i).calcCureRate(cure,T,Tg,cureRate);
    }
}

void Foam::kineticModelList::calcMaterialState
(
    const volScalarField& cureEff,
	const volScalarField& T,
	const volScalarField& TgEff,
	volScalarField& materialState
)
{
    forAll(*this, i)
    {
        this->operator[](i).calcMaterialState(cureEff,T,TgEff,materialState);
    }
}

void Foam::kineticModelList::calcMaterialState
(
    const volScalarField& cureEff,
	volScalarField& materialState
)
{
    forAll(*this, i)
    {
        this->operator[](i).calcMaterialState(cureEff,materialState);
    }
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const kineticModelList& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
