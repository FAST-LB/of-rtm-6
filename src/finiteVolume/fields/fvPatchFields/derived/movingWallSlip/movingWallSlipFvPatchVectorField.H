/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::movingWallSlipFvPatchVectorField

Group
    grpWallBoundaryConditions

Description


This *modified* boundary condition was created by Louis Gagnon ( louis.gagnon@polimi.it )
and is intented to be used as a moving wall boundary condition for AMI interface cases


    This boundary condition provides a velocity condition for cases with
    moving walls.  In addition, it should also be applied to 'moving' walls
    for moving reference frame (MRF) calculations.


    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        U            | velociy field name      | no          | U
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            movingWallSlip;
        U               U;
        value           uniform 0; // initial value
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchVectorField
    Foam::MRFZone

SourceFiles
    movingWallSlipFvPatchVectorField.C
    

\*---------------------------------------------------------------------------*/

#ifndef movingWallSlipFvPatchVectorField_H
#define movingWallSlipFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class movingWallSlipFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingWallSlipFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

       //- Name of velocity field
       word UName_;


public:

    //- Runtime type information
    TypeName("movingWallSlip");


    // Constructors

        //- Construct from patch and internal field
        movingWallSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingWallSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given movingWallSlipFvPatchVectorField
        //  onto a new patch
        movingWallSlipFvPatchVectorField
        (
            const movingWallSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingWallSlipFvPatchVectorField
        (
            const movingWallSlipFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new movingWallSlipFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingWallSlipFvPatchVectorField
        (
            const movingWallSlipFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new movingWallSlipFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
