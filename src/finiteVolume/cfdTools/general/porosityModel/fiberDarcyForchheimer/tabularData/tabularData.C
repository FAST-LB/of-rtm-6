/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "tabularData.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tabularData, 0);
    addToRunTimeSelectionTable(permeabilityLaw, tabularData, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::tabularData::tabularData
(
    const fvMesh& mesh,
    const word& name,
    const label& cellZoneID,
    const dictionary& dict
)
:
    permeabilityLaw(mesh, name, cellZoneID, dict),
    mesh_(mesh),
    cellZoneID_(cellZoneID),
    DA_("DA", dimensionSet(0, -2, 0, 0, 0), Zero),
    FO_("FO", dimensionSet(0, -1, 0, 0, 0), Zero),    
    messfvg_(0),
    messk1_(0),
    messk2_(0),
    messk3_(0),
    A1_(0),
    A2_(0),
    A3_(0),
    B1_(0),
    B2_(0),
    B3_(0),
    KK1(),
    KK2(),
    KK3(),
    interpolationType_("default")
{   
    const dictionary& tabularDataDict(dict.subDict("tabularDataCoeffs"));
    
    scalar t=0;
        if
        (
            !tabularDataDict.readIfPresent("interpolation", interpolationType_)
            || !tabularDataDict.readIfPresent("FVF", messfvg_)
            || !tabularDataDict.readIfPresent("K1", messk1_)
            || !tabularDataDict.readIfPresent("K2", messk2_)
            || !tabularDataDict.readIfPresent("K3", messk3_)
        )
        {
            FatalIOErrorIn
            (
                    "Foam::porosityModel::fiberDarcyForchheimer"
                    "(const keyType&, const fvMesh&, const dictionary&)",
                    dict_
            )   << "entries missing in the porosityProperties.fiberDarcyForchheimerCoeffs.measuredData dict. The following have to be defined:\n"
            << "- interpolation, "
            << "FVF, "
            << "K1, "
            << "K2, "
            << "K3 "
            << exit(FatalIOError);
        }  
        if 
        (
            messfvg_.size() != messk1_.size() 
            || messfvg_.size() != messk2_.size()
            || messfvg_.size() != messk3_.size()
        )
        {        
            FatalIOErrorIn
            (
                    "Foam::porosityModel::fiberDarcyForchheimer"
                    "(const keyType&, const fvMesh&, const dictionary&)",
                    dict_
            )   << "entries do not have the same number of values\n"
            << "FVF, "
            << "K1, "
            << "K2, "
            << "K3 "
            << exit(FatalIOError);    
        }
        if (interpolationType_ == word ("exponential"))
        {
            A1_=messfvg_;
            A2_=messfvg_;
            A3_=messfvg_;
            B1_=messfvg_;
            B2_=messfvg_;
            B3_=messfvg_;
            for(scalar m=0;m<messfvg_.size()-1; m++)
            {

                B1_[m]=log(messk1_[m+1]/messk1_[m])/(messfvg_[m+1]-messfvg_[m]);
                A1_[m]=messk1_[m]/exp(B1_[m]*messfvg_[m]);

                B2_[m]=log(messk2_[m+1]/messk2_[m])/(messfvg_[m+1]-messfvg_[m]);
                A2_[m]=messk2_[m]/exp(B2_[m]*messfvg_[m]);

                B3_[m]=log(messk3_[m+1]/messk3_[m])/(messfvg_[m+1]-messfvg_[m]);
                A3_[m]=messk3_[m]/exp(B3_[m]*messfvg_[m]);

                t=t+1;
            }

            messfvg_[0] = 0;
            messfvg_[t] = 1;
            A1_[t] = 0;
            A2_[t] = 0;
            A3_[t] = 0;
            B1_[t] = 0;
            B2_[t] = 0;
            B3_[t] = 0;

            DA_.value().xx()=1;
            DA_.value().yy()=1;
            DA_.value().zz()=1;
            FO_.value().xx()=1;
            FO_.value().yy()=1;
            FO_.value().zz()=1;
        }
        else if (interpolationType_ == word ("linear") )
        {
            A1_=messfvg_;
            A2_=messfvg_;
            A3_=messfvg_;
            B1_=messfvg_;
            B2_=messfvg_;
            B3_=messfvg_;
            for(scalar m=0;m<messfvg_.size()-1; m++)
            {

                B1_[m]=(messk1_[m]-messk1_[m+1])/(messfvg_[m]-messfvg_[m+1]);
                A1_[m]=messk1_[m]-(B1_[m]*messfvg_[m]);

                B2_[m]=(messk2_[m]-messk2_[m+1])/(messfvg_[m]-messfvg_[m+1]);
                A2_[m]=messk2_[m]-(B2_[m]*messfvg_[m]);

                B3_[m]=(messk3_[m]-messk3_[m+1])/(messfvg_[m]-messfvg_[m+1]);
                A3_[m]=messk3_[m]-(B3_[m]*messfvg_[m]);

                t=t+1;
            }

            messfvg_[0] = 0;
            messfvg_[t] = 1;
            A1_[t] = 0;
            A2_[t] = 0;
            A3_[t] = 0;
            B1_[t] = 0;
            B2_[t] = 0;
            B3_[t] = 0;

            DA_.value().xx()=1;
            DA_.value().yy()=1;
            DA_.value().zz()=1;
            FO_.value().xx()=1;
            FO_.value().yy()=1;
            FO_.value().zz()=1;
        }
        else if (interpolationType_ == word ("spline") )
        {
            std::vector<double> mfvg, mk1, mk2, mk3;
            forAll (messfvg_, m)
            {
                mfvg.push_back(messfvg_[m]);
                mk1.push_back(messk1_[m]);
                mk2.push_back(messk2_[m]);
                mk3.push_back(messk3_[m]);
            }
            KK1.set_points(mfvg,mk1);
            KK2.set_points(mfvg,mk2);
            KK3.set_points(mfvg,mk3);

            DA_.value().xx()=1;
            DA_.value().yy()=1;
            DA_.value().zz()=1;
            FO_.value().xx()=1;
            FO_.value().yy()=1;
            FO_.value().zz()=1;
        }
        else
        {
            FatalIOErrorIn
            (
                    "Foam::porosityModel::fiberDarcyForchheimer"
                    "(const keyType&, const fvMesh&, const dictionary&)",
                    dict_
            )   << "wrong entry in the porosityProperties dict. "<< interpolationType_ <<" is no valid interpolation type. Valid types are:\n"
            << "exponential, "
            << "linear, "
            << "spline "
            << exit(FatalIOError);
        }     
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tabularData::~tabularData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tensor Foam::tabularData::D(const scalar& FVF_, const label& j)
{    
    tensor D_ = Zero;

    if (interpolationType_ == word ("exponential"))
    {
        forAll (messfvg_, m)
        {
            
            if (FVF_>=messfvg_[m] && FVF_<messfvg_[m+1])
            {                
                D_.xx() = scalar(1) / (A1_[m]*exp(B1_[m]*FVF_));    
                D_.yy() = scalar(1) / (A2_[m]*exp(B2_[m]*FVF_));    
                D_.zz() = scalar(1) / (A3_[m]*exp(B3_[m]*FVF_));                    
            }                    
        }
    }
    else if (interpolationType_ == word ("linear"))
    {
        forAll (messfvg_, m)
        {
        
            if (FVF_>=messfvg_[m] && FVF_<messfvg_[m+1])
            {                
                D_.xx() = scalar(1) / (A1_[m]+B1_[m]*FVF_);    
                D_.yy() = scalar(1) / (A2_[m]+B2_[m]*FVF_);    
                D_.zz() = scalar(1) / (A3_[m]+B3_[m]*FVF_);
            }                    
        }            
    }
    else if (interpolationType_ == word ("spline"))
    {
        D_.xx() = scalar(1) / KK1(FVF_);
        D_.yy() = scalar(1) / KK2(FVF_);            
        D_.zz() = scalar(1) / KK3(FVF_);            
    }    
    
    return D_;
}

Foam::tmp<Foam::tensorField> Foam::tabularData::D(const scalarField& FVF_)
{      
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField D_(cells.size());

    forAll(cells, j)
    {
                
        D_[j] = Zero;

        if (interpolationType_ == word ("exponential"))
        {
            forAll (messfvg_, m)
            {
            
                if (FVF_[j]>=messfvg_[m] && FVF_[j]<messfvg_[m+1])
                {                
                    D_[j].xx() = scalar(1) / (A1_[m]*exp(B1_[m]*FVF_[j]));    
                    D_[j].yy() = scalar(1) / (A2_[m]*exp(B2_[m]*FVF_[j]));    
                    D_[j].zz() = scalar(1) / (A3_[m]*exp(B3_[m]*FVF_[j]));                    
                }                    
            }
        }
        else if (interpolationType_ == word ("linear"))
        {
            forAll (messfvg_, m)
            {
            
                if (FVF_[j]>=messfvg_[m] && FVF_[j]<messfvg_[m+1])
                {                
                    D_[j].xx() = scalar(1) / (A1_[m]+B1_[m]*FVF_[j]);    
                    D_[j].yy() = scalar(1) / (A2_[m]+B2_[m]*FVF_[j]);    
                    D_[j].zz() = scalar(1) / (A3_[m]+B3_[m]*FVF_[j]);
                }                    
            }            
        }
        else if (interpolationType_ == word ("spline"))
        {
            D_[j].xx() = scalar(1) / KK1(FVF_[j]);
            D_[j].yy() = scalar(1) / KK2(FVF_[j]);            
            D_[j].zz() = scalar(1) / KK3(FVF_[j]);            
        } 
    }
    
    return tmp<tensorField>(D_);
}

Foam::tmp<Foam::tensorField> Foam::tabularData::D()
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField D_(cells.size());
    forAll(cells, j)
    {
        D_[j] = Zero;
    }
    return tmp<tensorField>(D_);
}

Foam::tensor Foam::tabularData::F(const scalar& FVF_, const label& j)
{
    return tensor::zero;
}

Foam::tmp<Foam::tensorField> Foam::tabularData::F(const scalarField& FVF_)
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField F_(cells.size());
    forAll(cells, j)
    {
        F_[j] = Zero;
    }
    
    return tmp<tensorField>(F_);
}

Foam::tmp<Foam::tensorField> Foam::tabularData::F()
{
    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    
    tensorField F_(cells.size());
    forAll(cells, j)
    {
        F_[j] = Zero;
    }

    return tmp<tensorField>(F_);
}
// ************************************************************************* //
