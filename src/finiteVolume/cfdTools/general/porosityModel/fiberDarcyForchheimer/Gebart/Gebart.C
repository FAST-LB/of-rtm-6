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

#include "Gebart.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Gebart, 0);
    addToRunTimeSelectionTable(permeabilityLaw, Gebart, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::Gebart::Gebart
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
    vf0_(0.5),
    vfmax_(0.9),
    constDA_(0.0),
    constDB1_(0.0),
    constDB2_(0.0),
    constFA_(0.0),
    constFB1_(0.0), 
    constFB2_(0.0)
{    
    const dictionary& GebartDict(dict.subDict("GebartCoeffs"));
    if
    (
        !GebartDict.readIfPresent("vf0", vf0_)
        || !GebartDict.readIfPresent("vfmax", vfmax_)                
    )
    {
        FatalIOErrorIn
        (
            "Foam::porosityModel::fiberDarcyForchheimer"
            "(const keyType&, const fvMesh&, const dictionary&)",
            GebartDict
        )   << "model constants missing in the porosityProperties.fiberDarcyForchheimerCoeffs.Gebart dict. The following have to be defined:\n"
        << "- reference fiber volume fraction (0 - 1): vf0\n"
        << "- maximum fiber volume fraction (0 - 1): vfmax\n"            
        << exit(FatalIOError);
    }

    dimensionedVector d(vector::zero);
    if (GebartDict.readIfPresent("d", d))
    {        
        if (DA_.dimensions() != d.dimensions())
        {
            FatalIOErrorIn
            (
                "Foam::porosityModel::fiberDarcyForchheimer"
                "(const keyType&, const fvMesh&, const dictionary&)",
                dict_
            )   << "incorrect dimensions for d: " << d.dimensions()
            << " should be " << DA_.dimensions()
            << exit(FatalIOError);
        }

        DA_.value().xx() = d.value().x()+SMALL;
        DA_.value().yy() = d.value().y()+SMALL;
        DA_.value().zz() = d.value().z()+SMALL;
    }
    else 
    {
        FatalIOErrorIn
        (
            "Foam::porosityModel::fiberDarcyForchheimer"
            "(const keyType&, const fvMesh&, const dictionary&)",
            dict_
        )   << "Darcy vector d has to be defined. "
        << exit(FatalIOError);
    }   
        
    // constA = 8 * Rf^2 / c, see Gebart (1992)
    constDA_ = pow(vf0_, 2) / ( DA_.value().xx() * pow((1 - vf0_), 3) );

    // constB = C1 * Rf^2, see Gebart (1992)
    // constB1: constB in equation for in-plane transversal permeability
    // constB2: constB in equation for out-of-plane transversal permeability
    constDB1_ = scalar(1) / (DA_.value().yy() * pow( (sqrt(vfmax_/vf0_) - 1), scalar(2.5)) );
    constDB2_ = scalar(1) / (DA_.value().zz() * pow( (sqrt(vfmax_/vf0_) - 1), scalar(2.5)) );    

    dimensionedVector f(vector::zero);
    if (GebartDict.readIfPresent("f", f))
    {
        if (FO_.dimensions() != f.dimensions())
        {
            FatalIOErrorIn
            (
                "Foam::porosityModel::fiberDarcyForchheimere"
                "(const keyType&, const fvMesh&, const dictionary&)",
                GebartDict
            )   << "incorrect dimensions for f: " << f.dimensions()
            << " should be " << FO_.dimensions()
            << exit(FatalIOError);
        }

        // leading 0.5 is from 1/2 * rho
        FO_.value().xx() = 0.5*f.value().x()+SMALL;
        FO_.value().yy() = 0.5*f.value().y()+SMALL;
        FO_.value().zz() = 0.5*f.value().z()+SMALL;
    }
    else 
    {
        FatalIOErrorIn
        (
            "Foam::porosityModel::fiberDarcyForchheimer"
            "(const keyType&, const fvMesh&, const dictionary&)",
            dict_
        )   << "Forchheimer vector f has to be defined. "
        << exit(FatalIOError);
    }
    
    // constA = 8 * Rf^2 / c, see Gebart (1992)
    constFA_ = pow(vf0_, 2) / (FO_.value().xx() * pow((1 - vf0_), 3) );

    // constB = C1 * Rf^2, see Gebart (1992)
    // constB1: constB in equation for in-plane transversal permeability
    // constB2: constB in equation for out-of-plane transversal permeability
    constFB1_ = scalar(1) / (FO_.value().yy() * pow( (sqrt(vfmax_/vf0_) - 1), scalar(2.5)) );
    constFB2_ = scalar(1) / (FO_.value().zz() * pow( (sqrt(vfmax_/vf0_) - 1), scalar(2.5)) );

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Gebart::~Gebart()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tensor Foam::Gebart::D(const scalar& FVF_, const label& j)
{ 
    tensor D_(Zero);

    D_.xx() = scalar(1) / (constDA_ * pow((1 - FVF_), 3) / pow(FVF_, 2));
    D_.yy() = scalar(1) / (constDB1_ * pow( (sqrt(vfmax_/FVF_) - 1), scalar(2.5)));
    D_.zz() = scalar(1) / (constDB2_ * pow( (sqrt(vfmax_/FVF_) - 1), scalar(2.5)));                       
    
    return D_;
}

Foam::tmp<Foam::tensorField> Foam::Gebart::D(const scalarField& FVF_)
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField D_(cells.size());

    forAll(cells, j)
    {        
        D_[j] = Zero;
        D_[j].xx() = scalar(1) / (constDA_ * pow((1 - FVF_[j]), 3) / pow(FVF_[j], 2));
        D_[j].yy() = scalar(1) / (constDB1_ * pow( (sqrt(vfmax_/FVF_[j]) - 1), scalar(2.5)));
        D_[j].zz() = scalar(1) / (constDB2_ * pow( (sqrt(vfmax_/FVF_[j]) - 1), scalar(2.5)));                       
    }

    return tmp<tensorField>(D_);
}

Foam::tmp<Foam::tensorField> Foam::Gebart::D()
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField D_(cells.size());
    forAll(cells, j)
    {
        D_[j] = Zero;
    }
    return tmp<tensorField>(D_);
}

Foam::tensor Foam::Gebart::F(const scalar& FVF_, const label& j)
{
    tensor F_(Zero);
         
    F_.xx() = scalar(1) / (constFA_ * pow((1 - FVF_), 3) / pow(FVF_, 2));
    F_.yy() = scalar(1) / (constFB1_ * pow( (sqrt(vfmax_/FVF_) - 1), scalar(2.5)));
    F_.zz() = scalar(1) / (constFB2_ * pow( (sqrt(vfmax_/FVF_) - 1), scalar(2.5)));                       
    
    return F_;            
}

Foam::tmp<Foam::tensorField> Foam::Gebart::F(const scalarField& FVF_)
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField F_(cells.size());
    
    forAll(cells, j)
    {        
        F_[j] = Zero;
        F_[j].xx() = scalar(1) / (constFA_ * pow((1 - FVF_[j]), 3) / pow(FVF_[j], 2));
        F_[j].yy() = scalar(1) / (constFB1_ * pow( (sqrt(vfmax_/FVF_[j]) - 1), scalar(2.5)));
        F_[j].zz() = scalar(1) / (constFB2_ * pow( (sqrt(vfmax_/FVF_[j]) - 1), scalar(2.5)));                       
    }

    return tmp<tensorField>(F_);
}

Foam::tmp<Foam::tensorField> Foam::Gebart::F()
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField F_(cells.size());
    forAll(cells, j)
    {
        F_[j] = Zero;
    }
    Info << F_ << endl;
    return tmp<tensorField>(F_);
}
// ************************************************************************* //
