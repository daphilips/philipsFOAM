/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

template<class ParcelType>
inline
Foam::ReactingParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    pMin_(0.0),
    constantVolume_(false),
    Tvap_(0.0),
    Tbp_(0.0)
{}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    pMin_(cp.pMin_),
    constantVolume_(cp.constantVolume_),
    Tvap_(cp.Tvap_),
    Tbp_(cp.Tbp_)
{}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict,
    const bool readFields
)
:
    ParcelType::constantProperties(parentDict, readFields),
    pMin_(0.0),
    constantVolume_(false),
    Tvap_(0.0),
    Tbp_(0.0)
{
    if (readFields)
    {
        this->dict().lookup("pMin") >> pMin_;
        this->dict().lookup("constantVolume") >> constantVolume_;
        this->dict().lookup("Tvap") >> Tvap_;
        this->dict().lookup("Tbp") >> Tbp_;
    }
}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::constantProperties::constantProperties
(
    const label parcelTypeId,
    const scalar rhoMin,
    const scalar rho0,
    const scalar minParticleMass,
    const scalar youngsModulus,
    const scalar poissonsRatio,
    const scalar T0,
    const scalar TMin,
    const scalar Cp0,
    const scalar epsilon0,
    const scalar f0,
    const scalar Pr,
    const scalar pMin,
    const Switch& constantVolume,
    const scalar Tvap,
    const scalar Tbp
)
:
    ParcelType::constantProperties
    (
        parcelTypeId,
        rhoMin,
        rho0,
        minParticleMass,
        youngsModulus,
        poissonsRatio,
        T0,
        TMin,
        Cp0,
        epsilon0,
        f0,
        Pr
    ),
    pMin_(pMin),
    constantVolume_(constantVolume),
    Tvap_(Tvap),
    Tbp_(Tbp)
{}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(mesh, position, cellI, tetFaceI, tetPtI),
    mass0_(0.0),
    Y_(0),
    pc_(0.0)
{}


template<class ParcelType>
inline Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const scalarField& Y0,
    const constantProperties& constProps
)
:
    ParcelType
    (
        mesh,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        constProps
    ),
    mass0_(0.0),
    Y_(Y0),
    pc_(0.0)
{
    // Set initial parcel mass
    mass0_ = this->mass();
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ReactingParcel<ParcelType>::constantProperties::pMin() const
{
    return pMin_;
}


template<class ParcelType>
inline Foam::Switch
Foam::ReactingParcel<ParcelType>::constantProperties::constantVolume() const
{
    return constantVolume_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingParcel<ParcelType>::constantProperties::Tvap() const
{
    return Tvap_;
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingParcel<ParcelType>::constantProperties::Tbp() const
{
    return Tbp_;
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::mass0() const
{
    return mass0_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingParcel<ParcelType>::Y() const
{
    return Y_;
}


template<class ParcelType>
inline Foam::scalar Foam::ReactingParcel<ParcelType>::pc() const
{
    return pc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::pc()
{
    return pc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ReactingParcel<ParcelType>::mass0()
{
    return mass0_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingParcel<ParcelType>::Y()
{
    return Y_;
}


// ************************************************************************* //
