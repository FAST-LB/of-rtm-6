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
#include "fiberDarcyForchheimer.H"
#include "zeroGradientFvPatchFields.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(fiberDarcyForchheimer, 0);
        addToRunTimeSelectionTable(porosityModel, fiberDarcyForchheimer, mesh);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::fiberDarcyForchheimer::fiberDarcyForchheimer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName  
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),    
    K_(),
    permeabilityLawPtrs_(cellZoneIDs_.size()),
    FVF_(cellZoneIDs_.size()),
    D_(cellZoneIDs_.size()),
    F_(cellZoneIDs_.size()),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho")),
    muName_(coeffs_.lookupOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffs_.lookupOrDefault<word>("nu", "nu")),
    drapeDir(mesh_.thisDb().lookupObject<volVectorField>("drapeDirection").internalField()),
    thicknessDir(mesh_.thisDb().lookupObject<volVectorField>("drapeDirection").internalField())
{  
    //Initialize permeability models    
    forAll(cellZoneIDs_, zoneI)
    {
        permeabilityLawPtrs_= permeabilityLaw::New(mesh_, name, zoneI, coeffs_);    
    }    
    calcThicknessDirection();  
    calcTransformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::fiberDarcyForchheimer::~fiberDarcyForchheimer()
{    
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volTensorField> Foam::porosityModels::fiberDarcyForchheimer::K()
{
    volTensorField permeability(*K_);
    forAll(cellZoneIDs_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
        forAll(cells, j)
        {            
            if (D_[zoneI][j].xx() > 1)
                permeability.ref()[cells[j]].xx() = 1 / D_[zoneI][j].xx();
            if (D_[zoneI][j].xy() > 1)
                permeability.ref()[cells[j]].xy() = 1 / D_[zoneI][j].xy();
            if (D_[zoneI][j].xz() > 1)
                permeability.ref()[cells[j]].xz() = 1 / D_[zoneI][j].xz();
            if (D_[zoneI][j].yx() > 1)
                permeability.ref()[cells[j]].yx() = 1 / D_[zoneI][j].yx();
            if (D_[zoneI][j].yy() > 1)            
                permeability.ref()[cells[j]].yy() = 1 / D_[zoneI][j].yy();
            if (D_[zoneI][j].yz() > 1)
                permeability.ref()[cells[j]].yz() = 1 / D_[zoneI][j].yz();
            if (D_[zoneI][j].zx() > 1)
                permeability.ref()[cells[j]].zx() = 1 / D_[zoneI][j].zx();
            if (D_[zoneI][j].zy() > 1)
                permeability.ref()[cells[j]].zy() = 1 / D_[zoneI][j].zy();
            if (D_[zoneI][j].zz() > 1)
                permeability.ref()[cells[j]].zz() = 1 / D_[zoneI][j].zz();
        }
    }
    
    return tmp<volTensorField>
    (
        new volTensorField
        (            
            permeability
        )
    );
}

//Calculate thickness direction for permeability rotation
void Foam::porosityModels::fiberDarcyForchheimer::calcThicknessDirection()
{   

    bool calcThickness(false);
    coeffs_.readIfPresent("calcThickness", calcThickness);
    if(calcThickness == true)
    {
        Info << "    Calculating local thickness direction for permeability rotation." << endl;
        
        //read wall patches from porosityProperties dictionary
        label patchID = mesh_.boundaryMesh().findPatchID(coeffs_.lookup("thicknessWall"));
        if ( patchID < 0 )
        {    
            FatalIOErrorIn
            (
                    "Foam::porosityModel::fiberDarcyForchheimer"
                    "(const keyType&, const fvMesh&, const dictionary&)",
                    dict_
            )   << "wrong entry in the porosityProperties dict. thicknessWall entry is not valid: "
                <<  coeffs_.lookup("thicknessWall") << " does not exist :\n"
            << exit(FatalIOError);    
        }        
        
        const labelList& faces = mesh_.boundaryMesh()[patchID].faceCells();   
        
        forAll(cellZoneIDs_, zoneI)
        {    
            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
            List <coordinateSystem> localCOS(cells.size(), coordinateSystem());            

            forAll(cells, j)
            { 
                //search for minimum distance to patch 1
                vector minDistance = vector(1,1,1);
                int k = 0;
                forAll(faces, i)
                {
                    vector distance = mesh_.Cf().boundaryField()[patchID][i]-mesh_.C()[cells[j]];
                    if (mag(minDistance)>mag(distance))
                    {
                        minDistance = distance;
                        k = i;
                    }
                } 
                
                thicknessDir[cells[j]] = mesh_.Sf().boundaryField()[patchID][k];
                //normalize
                thicknessDir[cells[j]] /= mag(thicknessDir[cells[j]]);
                               
                //normalize
                drapeDir[cells[j]] /= mag(drapeDir[cells[j]]);               
            }        
        }       
    }
    else
    {
        //Add option to give vector in coeff dict
        Info << "    Using (0 0 1) as thickness direction for permeability rotation." << endl;
   
        forAll(cellZoneIDs_, zoneI)
        {    
            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];
          
            forAll(cells, j)
            {
                thicknessDir[cells[j]] = vector(0,0,1);
            }
        }
    }
}

//Transform permeability using fiber orientations     
void Foam::porosityModels::fiberDarcyForchheimer::calcTransformModelData()
{ 

    forAll(cellZoneIDs_, zoneI)
    { 
        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];   
        
        D_[zoneI].setSize(cells.size());
        F_[zoneI].setSize(cells.size());
        FVF_[zoneI].setSize(cells.size());           
                
        List<coordinateSystem> localCOS(cells.size(), coordinateSystem());       
        
        forAll(cells, j)
        {           
            FVF_[zoneI][j] = mesh_.thisDb().lookupObject<volScalarField>("fiberVolFraction").internalField()[cells[j]];
            
            localCOS[j] = coordinateSystem("localCOS", vector(0,0,0), thicknessDir[cells[j]], drapeDir[cells[j]]);            
            const coordinateRotation& R = localCOS[j].R();         
            
            D_[zoneI][j] = R.transformTensor(permeabilityLawPtrs_[zoneI]->D(FVF_[zoneI][j],j));
            F_[zoneI][j] = R.transformTensor(permeabilityLawPtrs_[zoneI]->F(FVF_[zoneI][j],j));           
        }
        
    }
    
    K_ = const_cast<volTensorField*>(&mesh_.thisDb().lookupObject<volTensorField>("permeability"));
    *K_ = K();
    K_->correctBoundaryConditions();
}    

void Foam::porosityModels::fiberDarcyForchheimer::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), Zero);
    const scalarField& V = mesh_.V();

    apply(Udiag, Usource, V, rho, mu, U);

    force = Udiag*U - Usource;
}

//explicit porosity
void Foam::porosityModels::fiberDarcyForchheimer::correct
(
    fvVectorMatrix& UEqn
) const
{
   
    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);

        if (mesh_.foundObject<volScalarField>(muName))
        {
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, rho, mu, U);
        }
        else
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, rho, rho*nu, U);
        }
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, geometricOneField(), nu, U);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, geometricOneField(), mu/rho, U);
        }
    }
}

//explicit porosity
void Foam::porosityModels::fiberDarcyForchheimer::correct
(
    fvVectorMatrix& UEqn, 
    const volScalarField& cellIbMask
) const
{
   
    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const scalarField& cellIbMask_ = cellIbMask.internalField();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);

        if (mesh_.foundObject<volScalarField>(muName))
        {
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, rho, mu, U, cellIbMask_);
        }
        else
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, rho, rho*nu, U, cellIbMask_);
        }
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, geometricOneField(), nu, U, cellIbMask_);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, geometricOneField(), mu/rho, U, cellIbMask_);
        }
    }
}

void Foam::porosityModels::fiberDarcyForchheimer::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& cellIbMask,
    const volVectorField& UFibers
) const
{
    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const scalarField& cellIbMask_ = cellIbMask.internalField();
    const vectorField& UFibers_ = UFibers.internalField();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);

        if (mesh_.foundObject<volScalarField>(muName))
        {
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, rho, mu, U, cellIbMask_, UFibers_);
        }
        else
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, rho, rho*nu, U, cellIbMask_, UFibers_);
        }
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, geometricOneField(), nu, U, cellIbMask_, UFibers_);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, geometricOneField(), mu/rho, U, cellIbMask_, UFibers_);
        }
    }
}

//explicit porosity
void Foam::porosityModels::fiberDarcyForchheimer::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{

    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    apply(Udiag, Usource, V, rho, mu, U);
}

//implicit porosity
void Foam::porosityModels::fiberDarcyForchheimer::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);
        const volScalarField& mu = mesh_.lookupObject<volScalarField>(muName);

        apply(AU, rho, mu, U);
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(AU, geometricOneField(), nu, U);            
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(AU, geometricOneField(), mu/rho, U);
        }
    }
}


bool Foam::porosityModels::fiberDarcyForchheimer::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
