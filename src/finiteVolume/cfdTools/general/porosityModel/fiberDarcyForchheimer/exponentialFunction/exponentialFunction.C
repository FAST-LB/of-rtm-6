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

\*---------------------------------------------------------------------------*/

#include "exponentialFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(exponentialFunction, 0);
    addToRunTimeSelectionTable(permeabilityLaw, exponentialFunction, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::exponentialFunction::exponentialFunction
(
    const fvMesh& mesh,
    const word& name,
    const label& cellZoneID,
    const dictionary& dict
)
:
    permeabilityLaw(mesh, name, cellZoneID, dict),
    mesh_(mesh),
    cellZoneID_(cellZoneID),
    DA_("DA", dimensionSet(0, -2, 0, 0, 0), Zero),
    FO_("FO", dimensionSet(0, -1, 0, 0, 0), Zero),    
    C01_(0),
    C02_(0),
    C03_(0),
    C11_(0),
    C12_(0),
    C13_(0)
{   
    const dictionary& exponentialFunctionDict(dict.subDict("exponentialFunctionCoeffs"));
    if
        (
            !exponentialFunctionDict.readIfPresent("C01", C01_)
            || !exponentialFunctionDict.readIfPresent("C02", C02_)
            || !exponentialFunctionDict.readIfPresent("C03", C03_)
            || !exponentialFunctionDict.readIfPresent("C11", C11_)
            || !exponentialFunctionDict.readIfPresent("C12", C12_)
            || !exponentialFunctionDict.readIfPresent("C13", C13_)
        )
        {
            FatalIOErrorIn
            (
                    "Foam::porosityModel::fiberDarcyForchheimer"
                    "(const keyType&, const fvMesh&, const dictionary&)",
                    dict_
            )   << "entries missing in the porosityProperties.fiberDarcyForchheimerCoeffs.exponentialFunction dict. The following have to be defined:\n"
            << "C01, "
            << "C02, "
            << "C03, "
            << "C11, "
            << "C12, "
            << "C13 \n"
            << " to calculate K1=C01*exp(C11*fiberVolFraction), K2=C02*exp(C12*fiberVolFraction), K1=C03*exp(C13*fiberVolFraction)"
            << exit(FatalIOError);
        }


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::exponentialFunction::~exponentialFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tensor Foam::exponentialFunction::D(const scalar& FVF_, const label& j)
{ 
    tensor D_(Zero);
         
    D_.xx() = scalar(1) / (C01_*exp(C11_*FVF_));    
    D_.yy() = scalar(1) / (C02_*exp(C12_*FVF_));    
    D_.zz() = scalar(1) / (C03_*exp(C13_*FVF_));                 
    
    return D_;
}

Foam::tmp<Foam::tensorField> Foam::exponentialFunction::D(const scalarField& FVF_)
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField D_(cells.size());

    forAll(cells, j)
    {
        D_[j] = Zero;
                
        D_[j].xx() = scalar(1) / (C01_*exp(C11_*FVF_[j]));    
        D_[j].yy() = scalar(1) / (C02_*exp(C12_*FVF_[j]));    
        D_[j].zz() = scalar(1) / (C03_*exp(C13_*FVF_[j]));                
    }
    
    return tmp<tensorField>(D_);
}

Foam::tmp<Foam::tensorField> Foam::exponentialFunction::D()
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField D_(cells.size());
    forAll(cells, j)
    {
        D_[j] = Zero;
    }
    return tmp<tensorField>(D_);
}

Foam::tensor Foam::exponentialFunction::F(const scalar& FVF_, const label& j)
{
    return tensor::zero;
}

Foam::tmp<Foam::tensorField> Foam::exponentialFunction::F(const scalarField& FVF_)
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField F_(cells.size());
    forAll(cells, j)
    {
        F_[j] = Zero;
    }
    
    return tmp<tensorField>(F_);
}

Foam::tmp<Foam::tensorField> Foam::exponentialFunction::F()
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField F_(cells.size());
    forAll(cells, j)
    {
        F_[j] = Zero;
    }

    return tmp<tensorField>(F_);
}
// ************************************************************************* //
