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
    Foam::kineticModel

Description
    Top level model for kinetic models

    Alexander Bernath, KIT-FAST-LBT, June 2017
    alexander.bernath@kit.edu

SourceFiles
    kineticModel.C
    kineticModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef kineticModel_H
#define kineticModel_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "dimensionedVector.H"
#include "keyType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class kineticModel Declaration
\*---------------------------------------------------------------------------*/

class kineticModel
:
    public regIOobject
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
		kineticModel(const kineticModel&);

        //- Disallow default bitwise assignment
        void operator=(const kineticModel&);


protected:

    // Protected data

        //- Kinetic model name
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Dictionary used for model construction
        const dictionary dict_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Model coeffient cure degree at gelation (user-supplied) [-]
        dimensionedScalar cureGel_;

        //- Kinetic active flag
        bool active_;

        //- Name(s) of cell-zone
        keyType zoneName_;

        //- Cell zone IDs
        labelList cellZoneIDs_;

    // Protected Member Functions

//        //- Return label index
//        label fieldIndex(const label index) const;


public:

    //- Runtime type information
    TypeName("kineticModel");

    //- Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        kineticModel,
        mesh,
        (
            const word& modelName,
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cellZoneName
        ),
        (modelName, name, mesh, dict, cellZoneName)
    );

    //- Constructor
    kineticModel
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName = word::null
    );

    //- Return pointer to new kineticModel object created on the freestore
    //  from an Istream
    class iNew
    {
        //- Reference to the mesh database
        const fvMesh& mesh_;
        const word& name_;

    public:

        iNew
        (
            const fvMesh& mesh,
            const word& name
        )
        :
            mesh_(mesh),
            name_(name)
        {}

        autoPtr<kineticModel> operator()(Istream& is) const
        {
            const dictionary dict(is);

            return autoPtr<kineticModel>
            (
				kineticModel::New
                (
                    name_,
                    mesh_,
                    dict
                )
            );
        }
    };

    //- Selector
    static autoPtr<kineticModel> New
    (
        const word& name,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName = word::null
    );

    //- Destructor
    virtual ~kineticModel();


    // Member Functions

        //- Return const access to the kinetic model name
        inline const word& name() const;

        //- Return const access to the kinetic active flag
        inline bool active() const;

        //- Return const access to the cell zone IDs
        inline const labelList& cellZoneIDs() const;

        //- Calculate current cure degree from time step and cure rate (given by model)
        void calcCure
		(
			volScalarField& cure,
			const volScalarField& cureRate
		);

        //- Calculate material state
        void calcMaterialState
        (
			const volScalarField& cureEff,
            const volScalarField& T,
            const volScalarField& TgEff,
			volScalarField& materialState
        );

        //- Calculate material state (no Tg data available)
        void calcMaterialState
        (
			const volScalarField& cureEff,
			volScalarField& materialState
        );

        //- Calculate the cure rate
        virtual void calcCureRate
        (
			const volScalarField& cure,
            const volScalarField& T,
            const volScalarField& Tg,
			volScalarField& cureRate
        ) const = 0;

    // I-O

        //- Write
        virtual bool writeData(Ostream& os) const;

        //- Inherit read from regIOobject
        using regIOobject::read;

        //- Read kinetic dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "kineticModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
