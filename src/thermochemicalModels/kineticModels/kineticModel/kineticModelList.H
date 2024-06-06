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
    Foam::kineticModelList

Description
    List container for kinetic models

SourceFiles
    kineticModelList.C

\*---------------------------------------------------------------------------*/

#ifndef kineticModelList_H
#define kineticModelList_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "kineticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class kineticModelList;
Ostream& operator<<(Ostream& os, const kineticModelList& models);

/*---------------------------------------------------------------------------*\
                      Class kineticModelList Declaration
\*---------------------------------------------------------------------------*/

class kineticModelList
:
    public PtrList<kineticModel>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        kineticModelList(const kineticModelList&);

        //- Disallow default bitwise assignment
        void operator=(const kineticModelList&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;


public:

    //- Constructor
    kineticModelList(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    ~kineticModelList();


    // Member Functions

        //- Return active status
        bool active(const bool active = false) const;

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Calculate cure degree
        void calcCure
        (
            volScalarField& cure,
            const volScalarField& cureRate
        );

        //- Calculate cure rate
        void calcCureRate
        (
            const volScalarField& cure,
            const volScalarField& T,
            const volScalarField& Tg,
			volScalarField& cureRate
        );

        //- Calculate material state
        void calcMaterialState
        (
            const volScalarField& cureEff,
            const volScalarField& T,
            const volScalarField& TgEff,
            volScalarField& materialState
        );

        //- Calculate material state
        void calcMaterialState
        (
            const volScalarField& cureEff,
            volScalarField& materialState
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
                const kineticModelList& models
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
