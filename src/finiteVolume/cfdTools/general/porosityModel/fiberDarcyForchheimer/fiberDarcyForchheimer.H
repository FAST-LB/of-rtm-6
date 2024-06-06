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
    Foam::fiberDarcyForchheimer

Description
    fiberDarcy-Forchheimer law porosity model, given by:

        \f[
            S = - (\mu d + \frac{\rho |U|}{2} f) U
        \f]

    where
    \vartable
        d        | Darcy coefficient [1/m2]
        f        | Forchheimer coefficient [1/m]
    \endvartable

    Since negative Darcy/Forchheimer parameters are invalid, they can be used
    to specify a multiplier (of the max component).
    
    Permeability K = 1/d
    Permeability is calculated with fiber volume fraction and fiber orientation and interpolation
    based on measuredData
       interpolation type exponential/linear/spline
       FVF( -values- )
       K1( -values- )
       K2( -values- )
       K3( -values- )
    based on Gebart
        vf0
        vfmax
        
    
    Fiber volume fraction
        on element defined by volumeScalarField fiberVolFraction
    or
        on porous zone by value fiberVolFraction
    
    Orientation
        on element defined by volVectorField drapeDirection 
    or
        on porous zone defined with the same notation as
        a co-ordinate system, but only a Cartesian co-ordinate system is valid
        


SourceFiles
    fiberDarcyForchheimer.C
    fiberDarcyForchheimerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fiberDarcyForchheimer_H
#define fiberDarcyForchheimer_H

#include "porosityModel.H"
#include "permeabilityLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                      Class fiberDarcyForchheimer Declaration
\*---------------------------------------------------------------------------*/

class fiberDarcyForchheimer
:
    public porosityModel
{
    // Private data
        
        //- Reference to stress field
        volTensorField* K_;

        //- Rheology law
        List<autoPtr<permeabilityLaw>> permeabilityLawPtrs_;
        
        List<scalarField> FVF_;

        //- Darcy coefficient - converted from dXYZ [1/m2]
        List<tensorField> D_;

        //- Forchheimer coefficient - converted from fXYZ [1/m]
        List<tensorField> F_;
        
        //- Name of density field
        word rhoName_;

        //- Name of dynamic viscosity field
        word muName_;

        //- Name of kinematic viscosity field
        word nuName_;
        
        //local values
        vectorField drapeDir;
        vectorField thicknessDir;



    // Private Member Functions

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;
        
        //- Apply
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U,
            const scalarField& cellIbMask
        ) const;
        
        //- Apply
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U,
            const scalarField& cellIbMask,
            const vectorField& UFibers
        ) const;

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            tensorField& AU,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        fiberDarcyForchheimer(const fiberDarcyForchheimer&);

        //- Disallow default bitwise assignment
        void operator=(const fiberDarcyForchheimer&);


public:

    //- Runtime type information
    TypeName("fiberDarcyForchheimer");

    //- Constructor
    fiberDarcyForchheimer
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~fiberDarcyForchheimer();


    // Member Functions

	//- Return permeability
        virtual tmp<volTensorField> K();
        
        virtual void calcThicknessDirection();

        //- Transform the model data wrt mesh changes
        virtual void calcTransformModelData();

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const;

        //- Add resistance
        virtual void correct(fvVectorMatrix& UEqn) const;
        
        //- Add resistance
        virtual void correct(fvVectorMatrix& UEqn, const volScalarField& cellIbMask) const;
        
        //- Add resistance
        virtual void correct(fvVectorMatrix& UEqn, const volScalarField& cellIbMask, const volVectorField& UFibers) const;

        //- Add resistance
        virtual void correct
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        //- Add resistance
        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;


    // I-O

        //- Write
        bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fiberDarcyForchheimerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
