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

Class
    Foam::TgModelList

Description
    List container for Tg models

SourceFiles
    TgModelList.C

\*---------------------------------------------------------------------------*/

#ifndef TgModelList_H
#define TgModelList_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "TgModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class TgModelList;
Ostream& operator<<(Ostream& os, const TgModelList& models);

/*---------------------------------------------------------------------------*\
                      Class TgModelList Declaration
\*---------------------------------------------------------------------------*/

class TgModelList
:
    public PtrList<TgModel>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        TgModelList(const TgModelList&);

        //- Disallow default bitwise assignment
        void operator=(const TgModelList&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;


public:

    //- Constructor
    TgModelList(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    ~TgModelList();


    // Member Functions

        //- Return active status
        bool active(const bool active = false) const;

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Calculate current glass transition temperature
        void calcTg
        (
            const volScalarField& cure,
            const volScalarField& T,
            volScalarField& Tg
        );

        // I-O

            //- Read dictionary
            bool read(const dictionary& dict);

            //- Write data to Ostream
            bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const TgModelList& models
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
