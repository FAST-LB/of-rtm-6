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

#include "HeleShaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HeleShaw, 0);
    addToRunTimeSelectionTable(permeabilityLaw, HeleShaw, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::HeleShaw::HeleShaw
(
    const fvMesh& mesh,
    const word& name,
    const label& cellZoneID,
    const dictionary& dict
)
:
    permeabilityLaw(mesh, name, cellZoneID, dict),
    mesh_(mesh),
    HeleShawDict(dict.subDict("HeleShawCoeffs")),
    cellZoneID_(cellZoneID),  
    H0_(HeleShawDict.lookupOrDefault("h0", 1.0)),
    Hdot_(HeleShawDict.lookupOrDefault("hdot", 0.0)),
    mu_(0.0),
    localGapHeight_(false),
    patchNames_(),
    faceZoneNames_(),
    patchIDs_(2),    
    faces1_(),
    faces2_(),
    height_(),
    timeValue_(0.0)
{    
    
    HeleShawDict.readIfPresent("localGapHeight", localGapHeight_);
    
    if
    (
        !HeleShawDict.readIfPresent("mu", mu_)        
    )
    {
        FatalIOErrorIn
        (
            "Foam::porosityModel::fiberDarcyForchheimer"
            "(const keyType&, const fvMesh&, const dictionary&)",
            HeleShawDict
        )   << "model constants missing in the porosityProperties.fiberDarcyForchheimerCoeffs.HeleShaw dict. The following have to be defined:\n"
        << "- viscosity scale factor: mu\n"
        << "- initial height and velocity: h0, hdot\n" 
        << "optional for local gap height:\n"
        << "- face set names: localGapHeight=true, patchNames\n"        
        << exit(FatalIOError);       
    }    
    
    if
    (
        HeleShawDict.readIfPresent("localGapHeight", localGapHeight_)
        && localGapHeight_ == true
    )
    {
        int nPatches = 0;
        
        //read patchNames from porosityProperties dictionary
        HeleShawDict.readIfPresent("patchNames", patchNames_);
        
        forAll(patchNames_,patchID)
        {        
            patchIDs_[patchID] = mesh_.boundaryMesh().findPatchID(patchNames_[patchID]);
            nPatches += 1;
        }
        if(nPatches < 2)
        {
            //read faceZoneNames from porosityProperties dictionary if present  
            HeleShawDict.readIfPresent("faceZoneNames", faceZoneNames_);
            forAll(faceZoneNames_,zoneID)
            {        
                patchIDs_[nPatches+zoneID] = mesh_.faceZones().findIndices(faceZoneNames_[zoneID])[zoneID];
                nPatches += 1;
            }            
        } 
        
        if(nPatches !=2)
        {
            FatalIOErrorIn
            (
                "Foam::porosityModel::fiberDarcyForchheimer"
                "(const keyType&, const fvMesh&, const dictionary&)",
                HeleShawDict
            )   << "Exactly 2 patches/faceZones have to be defined."
            << exit(FatalIOError);       
        }    
        
        faces1_ = mesh_.boundaryMesh()[patchIDs_[0]].faceCells();
        
        if(patchNames_.size() < 2)
        {
            faces2_ = mesh_.boundaryMesh()[patchIDs_[1]].faceCells();
        }
        else
        {
            faces2_ = mesh_.boundaryMesh()[patchIDs_[1]].faceCells();
        }
        
        height_ = H(mesh_.time().value());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HeleShaw::~HeleShaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//calculate height once each timestep
Foam::scalarField Foam::HeleShaw::H(const scalar& t)
{       
    
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    vectorField distances1(cells.size());
    vectorField distances2(cells.size());
    
            
        forAll(faces1_, k)
        {
            forAll(cells,j)
            {
                distances1[j] = mesh_.Cf().boundaryField()[patchIDs_[0]][k] - mesh_.C()[j];               
            }
        }
        forAll(faces2_, l)
        {
            forAll(cells,j)
            {
                distances2[j] = mesh_.Cf().boundaryField()[patchIDs_[1]][l] - mesh_.C()[j];
            }
        }       
    
    timeValue_ = t;
    
    const scalarField htemp = mag(distances1 - distances2);
    
    Info << "    New height: min = "<< min(htemp) <<" max = " << max(htemp) << endl;
        
    return htemp;
}

Foam::scalar Foam::HeleShaw::H(const scalar& t, const label& j)
{
    if(t != timeValue_)
    {
        height_ = H(mesh_.time().value());        
    }   
         
    return height_[j];
}

Foam::scalar Foam::HeleShaw::H(const label& j)
{  
    
    scalar gapHeight = 1.0;    
    
    
    if(localGapHeight_ == true)
    {  
        gapHeight = H(mesh_.time().value(), j);
    } 
    else
    {
        //Add option to give initial height and hdot in coeff dict         
        gapHeight = H0_ + mesh_.time().value()*Hdot_;               
        
    }    
    
    return gapHeight;

}

Foam::tensor Foam::HeleShaw::D(const scalar& FVF_, const label& j)
{    
    tensor D_(Zero);    
    const scalar H_ = H(j);
    
    D_.xx() = 12*mu_ / pow(H_, 2);
    D_.yy() = 12*mu_ / pow(H_, 2);
    D_.zz() = 12*mu_ / pow(H_, 2);                       
    
    //Info << D_ << endl;
    
    return D_;
}

Foam::tmp<Foam::tensorField> Foam::HeleShaw::D(const scalarField& FVF_)
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    tensorField D_(cells.size());

    forAll(cells, j)
    {      
        
        D_[j] = Zero;
        //D_[j].xx() = 12*mu_ / pow(H()[j], 2);
        //D_[j].yy() = 12*mu_ / pow(H()[j], 2);
        //D_[j].zz() = 12*mu_ / pow(H()[j], 2);                       
    }

    return tmp<tensorField>(D_);
}

Foam::tmp<Foam::tensorField> Foam::HeleShaw::D()
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField D_(cells.size());
    forAll(cells, j)
    {
        D_[j] = Zero;
    }
    return tmp<tensorField>(D_);
}

Foam::tensor Foam::HeleShaw::F(const scalar& FVF_, const label& j)
{
        
    return tensor::zero;            
}

Foam::tmp<Foam::tensorField> Foam::HeleShaw::F(const scalarField& FVF_)
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField F_(cells.size());
    
    forAll(cells, j)
    {        
        F_[j] = Zero;                            
    }

    return tmp<tensorField>(F_);
}

Foam::tmp<Foam::tensorField> Foam::HeleShaw::F()
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
