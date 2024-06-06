/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::CastroMacoskoTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    CastroMacoskoTransportI.H
    CastroMacoskoTransport.C

\*---------------------------------------------------------------------------*/

#ifndef CastroMacoskoTransport_H
#define CastroMacoskoTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class CastroMacoskoTransport;

template<class Thermo>
inline CastroMacoskoTransport<Thermo> operator+
(
    const CastroMacoskoTransport<Thermo>&,
    const CastroMacoskoTransport<Thermo>&
);

template<class Thermo>
inline CastroMacoskoTransport<Thermo> operator-
(
    const CastroMacoskoTransport<Thermo>&,
    const CastroMacoskoTransport<Thermo>&
);

template<class Thermo>
inline CastroMacoskoTransport<Thermo> operator*
(
    const scalar,
    const CastroMacoskoTransport<Thermo>&
);

template<class Thermo>
inline CastroMacoskoTransport<Thermo> operator==
(
    const CastroMacoskoTransport<Thermo>&,
    const CastroMacoskoTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const CastroMacoskoTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class CastroMacoskoTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class CastroMacoskoTransport
:
    public Thermo
{
    // Private data

        //- Castro-Macosko coefficients
        scalar B_, Tb_, C1_, C2_, cureGel_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline CastroMacoskoTransport
        (
            const Thermo& t,
            const scalar B,
            const scalar Tb,
            const scalar C1,
            const scalar C2,
            const scalar cureGel,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline CastroMacoskoTransport(const word&, const CastroMacoskoTransport&);

        //- Construct from Istream
        CastroMacoskoTransport(Istream&);

        //- Construct from dictionary
        CastroMacoskoTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<CastroMacoskoTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<CastroMacoskoTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<CastroMacoskoTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "CastroMacosko<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T, const scalar cure) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T, const scalar cure) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T, const scalar cure) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const CastroMacoskoTransport&);

        inline void operator+=(const CastroMacoskoTransport&);

        inline void operator-=(const CastroMacoskoTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend CastroMacoskoTransport operator+ <Thermo>
        (
            const CastroMacoskoTransport&,
            const CastroMacoskoTransport&
        );

        friend CastroMacoskoTransport operator- <Thermo>
        (
            const CastroMacoskoTransport&,
            const CastroMacoskoTransport&
        );

        friend CastroMacoskoTransport operator* <Thermo>
        (
            const scalar,
            const CastroMacoskoTransport&
        );

        friend CastroMacoskoTransport operator== <Thermo>
        (
            const CastroMacoskoTransport&,
            const CastroMacoskoTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const CastroMacoskoTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CastroMacoskoTransportI.H"

#ifdef NoRepository
    #include "CastroMacoskoTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
