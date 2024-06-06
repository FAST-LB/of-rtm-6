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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::CastroMacoskoTransport<Thermo>::CastroMacoskoTransport
(
    const Thermo& t,
    const scalar B,
    const scalar Tb,
    const scalar C1,
    const scalar C2,
    const scalar cureGel,
    const scalar Pr
)
:
    Thermo(t),
    B_(B),
    Tb_(Tb),
    C1_(C1),
    C2_(C2),
    cureGel_(cureGel),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::CastroMacoskoTransport<Thermo>::CastroMacoskoTransport
(
    const word& name,
    const CastroMacoskoTransport& cmt
)
:
    Thermo(name, cmt),
    B_(cmt.B_),
    Tb_(cmt.Tb_),
    C1_(cmt.C1_),
    C2_(cmt.C2_),
    cureGel_(cmt.cureGel_),
    rPr_(cmt.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::CastroMacoskoTransport<Thermo>>
Foam::CastroMacoskoTransport<Thermo>::clone() const
{
    return autoPtr<CastroMacoskoTransport<Thermo>>
    (
        new CastroMacoskoTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::CastroMacoskoTransport<Thermo>>
Foam::CastroMacoskoTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<CastroMacoskoTransport<Thermo>>
    (
        new CastroMacoskoTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::CastroMacoskoTransport<Thermo>>
Foam::CastroMacoskoTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<CastroMacoskoTransport<Thermo>>
    (
        new CastroMacoskoTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::CastroMacoskoTransport<Thermo>::mu
(
    const scalar p,
    const scalar T,
	const scalar cure
) const
{
    return B_*exp(Tb_/T)*pow(cureGel_/max(0.01,(cureGel_-max(0.0,cure))),(C1_+C2_*max(0.0,cure)));
//    return B_*Tb_/T*pow(cureGel_/(cureGel_-max(0.0,cure)),(C1_+C2_*max(0.0,cure)));
//    return scalar(1.0);//exp(Tb_/T);
}


template<class Thermo>
inline Foam::scalar Foam::CastroMacoskoTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T,
    const scalar cure
) const
{
    return this->Cp(p, T)*mu(p, T, cure)*rPr_;
//    return this->Cp(p, T)*scalar(0.1)*rPr_;
}


template<class Thermo>
inline Foam::scalar Foam::CastroMacoskoTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T,
    const scalar cure
) const
{
    return mu(p, T, cure)*rPr_;
//    return scalar(0.1)*rPr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::CastroMacoskoTransport<Thermo>::operator=
(
    const CastroMacoskoTransport<Thermo>& cmt
)
{
    Thermo::operator=(cmt);

    B_ = cmt.B_;
    Tb_ = cmt.Tb_;
    C1_ = cmt.C1_;
    C2_ = cmt.C2_;
    cureGel_ = cmt.cureGel_;
    rPr_ = cmt.rPr_;
}


template<class Thermo>
inline void Foam::CastroMacoskoTransport<Thermo>::operator+=
(
    const CastroMacoskoTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        B_ = molr1*B_ + molr2*st.B_;
        Tb_ = molr1*Tb_ + molr2*st.Tb_;
        C1_ = molr1*C1_ + molr2*st.C1_;
        C2_ = molr1*C2_ + molr2*st.C2_;
        cureGel_ = molr1*cureGel_ + molr2*st.cureGel_;
        rPr_ = 1.0/(molr1/rPr_ + molr2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::CastroMacoskoTransport<Thermo>::operator-=
(
    const CastroMacoskoTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        B_ = molr1*B_ - molr2*st.B_;
        Tb_ = molr1*Tb_ - molr2*st.Tb_;
        C1_ = molr1*C1_ - molr2*st.C1_;
        C2_ = molr1*C2_ - molr2*st.C2_;
        cureGel_ = molr1*cureGel_ - molr2*st.cureGel_;
        rPr_ = 1.0/(molr1/rPr_ - molr2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::CastroMacoskoTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::CastroMacoskoTransport<Thermo> Foam::operator+
(
    const CastroMacoskoTransport<Thermo>& cmt1,
    const CastroMacoskoTransport<Thermo>& cmt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(cmt1) + static_cast<const Thermo&>(cmt2)
    );

    if (mag(cmt1.nMoles()) + mag(cmt2.nMoles()) < SMALL)
    {
        return CastroMacoskoTransport<Thermo>
        (
            t,
            0,
            cmt1.rPr_
        );
    }
    else
    {
        scalar molr1 = cmt1.nMoles()/t.nMoles();
        scalar molr2 = cmt2.nMoles()/t.nMoles();

        return CastroMacoskoTransport<Thermo>
        (
            t,
            molr1*cmt1.B_ + molr2*cmt2.B_,
            molr1*cmt1.Tb_ + molr2*cmt2.Tb_,
            molr1*cmt1.C1_ + molr2*cmt2.C1_,
            molr1*cmt1.C2_ + molr2*cmt2.C2_,
            molr1*cmt1.cureGel_ + molr2*cmt2.cureGel_,
            1.0/(molr1/cmt1.rPr_ + molr2/cmt2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::CastroMacoskoTransport<Thermo> Foam::operator-
(
    const CastroMacoskoTransport<Thermo>& cmt1,
    const CastroMacoskoTransport<Thermo>& cmt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(cmt1) - static_cast<const Thermo&>(cmt2)
    );

    if (mag(cmt1.nMoles()) + mag(cmt2.nMoles()) < SMALL)
    {
        return CastroMacoskoTransport<Thermo>
        (
            t,
            0,
            cmt1.rPr_
        );
    }
    else
    {
        scalar molr1 = cmt1.nMoles()/t.nMoles();
        scalar molr2 = cmt2.nMoles()/t.nMoles();

        return CastroMacoskoTransport<Thermo>
        (
            t,
            molr1*cmt1.B_ - molr2*cmt2.B_,
            molr1*cmt1.Tb_ - molr2*cmt2.Tb_,
            molr1*cmt1.C1_ - molr2*cmt2.C1_,
            molr1*cmt1.C2_ - molr2*cmt2.C2_,
            molr1*cmt1.cureGel_ - molr2*cmt2.cureGel_,
            1.0/(molr1/cmt1.rPr_ - molr2/cmt2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::CastroMacoskoTransport<Thermo> Foam::operator*
(
    const scalar s,
    const CastroMacoskoTransport<Thermo>& cmt
)
{
    return CastroMacoskoTransport<Thermo>
    (
        s*static_cast<const Thermo&>(cmt),
        cmt.B_,
        cmt.Tb_,
        cmt.C1_,
        cmt.C2_,
        cmt.cureGel_,
        1.0/cmt.rPr_
    );
}


template<class Thermo>
inline Foam::CastroMacoskoTransport<Thermo> Foam::operator==
(
    const CastroMacoskoTransport<Thermo>& cmt1,
    const CastroMacoskoTransport<Thermo>& cmt2
)
{
    return cmt2 - cmt1;
}


// ************************************************************************* //
