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
    Foam::TgModel

Description
    Top level model for Tg models

	Alexander Bernath, KIT-FAST-LBT, June 2017
    alexander.bernath@kit.edu

SourceFiles
    TgModel.C
    TgModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef TgModel_H
#define TgModel_H

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
                        Class TgModel Declaration
\*---------------------------------------------------------------------------*/

class TgModel
:
    public regIOobject
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
		TgModel(const TgModel&);

        //- Disallow default bitwise assignment
        void operator=(const TgModel&);


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
    TypeName("TgModel");

    //- Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        TgModel,
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
    TgModel
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName = word::null
    );

    //- Return pointer to new TgModel object created on the freestore
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

        autoPtr<TgModel> operator()(Istream& is) const
        {
            const dictionary dict(is);

            return autoPtr<TgModel>
            (
				TgModel::New
                (
                    name_,
                    mesh_,
                    dict
                )
            );
        }
    };

    //- Selector
    static autoPtr<TgModel> New
    (
        const word& name,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName = word::null
    );

    //- Destructor
    virtual ~TgModel();


    // Member Functions

        //- Return const access to the Tg model name
        inline const word& name() const;

        //- Return const access to the Tg active flag
        inline bool active() const;

        //- Return const access to the cell zone IDs
        inline const labelList& cellZoneIDs() const;

        //- Calculate the current glass transition temperature
        virtual void calcTg
        (
			const volScalarField& cure,
            const volScalarField& T,
			volScalarField& Tg
        ) const = 0;

    // I-O

        //- Write
        virtual bool writeData(Ostream& os) const;

        //- Inherit read from regIOobject
        using regIOobject::read;

        //- Read Tg dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TgModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
