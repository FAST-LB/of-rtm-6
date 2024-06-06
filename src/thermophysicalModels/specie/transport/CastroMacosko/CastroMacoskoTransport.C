/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "CastroMacoskoTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::CastroMacoskoTransport<Thermo>::CastroMacoskoTransport(Istream& is)
:
    Thermo(is),
    B_(readScalar(is)),
    Tb_(readScalar(is)),
    C1_(readScalar(is)),
    C2_(readScalar(is)),
    cureGel_(readScalar(is)),
    rPr_(1.0/readScalar(is))
{
    is.check("CastroMacoskoTransport::CastroMacoskoTransport(Istream& is)");
}


template<class Thermo>
Foam::CastroMacoskoTransport<Thermo>::CastroMacoskoTransport(const dictionary& dict)
:
    Thermo(dict),
    B_(readScalar(dict.subDict("transport").lookup("B"))),
    Tb_(readScalar(dict.subDict("transport").lookup("Tb"))),
    C1_(readScalar(dict.subDict("transport").lookup("C1"))),
    C2_(readScalar(dict.subDict("transport").lookup("C2"))),
    cureGel_(readScalar(dict.subDict("transport").lookup("cureGel"))),
    rPr_(1.0/readScalar(dict.subDict("transport").lookup("Pr")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::CastroMacoskoTransport<Thermo>::CastroMacoskoTransport::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("B", B_);
    dict.add("Tb", Tb_);
    dict.add("C1", C1_);
    dict.add("C2", C2_);
    dict.add("cureGel", cureGel_);
    dict.add("Pr", 1.0/rPr_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const CastroMacoskoTransport<Thermo>& cmt)
{
    operator<<(os, static_cast<const Thermo&>(cmt));
    os << tab << cmt.B_ << tab << cmt.Tb_
	   << tab << cmt.C1_ << tab << cmt.C2_
	   << tab << cmt.cureGel_ << tab << 1.0/cmt.rPr_;

    os.check("Ostream& operator<<(Ostream&, const CastroMacoskoTransport&)");

    return os;
}


// ************************************************************************* //
