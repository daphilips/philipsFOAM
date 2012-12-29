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

#include "GGCInjection.H"
#include "mathematicalConstants.H"
#include "PackedBoolList.H"
#include "Switch.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::GGCInjection<CloudType>::GGCInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    positionsFile_(this->coeffDict().lookup("positionsFile")),
    positions_
    (
        IOobject
        (
            positionsFile_,
            owner.db().time().constant(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    diameters_(positions_.size()),
    injectorCells_(positions_.size(), -1),
    injectorTetFaces_(positions_.size(), -1),
    injectorTetPts_(positions_.size(), -1),
    U0_(this->coeffDict().lookup("U0")),
    min_d_(readScalar(this->coeffDict().lookup("min_diameter"))),
    max_d_(readScalar(this->coeffDict().lookup("max_diameter"))),
    min_rho_(readScalar(this->coeffDict().lookup("min_density"))),
    max_rho_(readScalar(this->coeffDict().lookup("max_density"))),
    particleRhos_(positions_.size())
{
    Switch ignoreOutOfBounds
    (
        this->coeffDict().lookupOrDefault("ignoreOutOfBounds", false)
    );

    label nRejected = 0;

    PackedBoolList keep(positions_.size(), true);

    forAll(positions_, pI)
    {
        if
        (
            !this->findCellAtPosition
            (
                injectorCells_[pI],
                injectorTetFaces_[pI],
                injectorTetPts_[pI],
                positions_[pI],
                !ignoreOutOfBounds
            )
        )
        {
            keep[pI] = false;

            nRejected++;
        }
    }

    if (nRejected > 0)
    {
        inplaceSubset(keep, positions_);
        inplaceSubset(keep, diameters_);
        inplaceSubset(keep, injectorCells_);
        inplaceSubset(keep, injectorTetFaces_);
        inplaceSubset(keep, injectorTetPts_);
        inplaceSubset(keep, particleRhos_);

        Info<< "    " << nRejected
            << " particles ignored, out of bounds." << endl;
    }

    // Construct parcel diameters
    scalar ddiam, drho;
    //scalar max_d_ = 1.0;
    //scalar min_d_ = 0.1;
    if (positions_.size() > 1){
      ddiam = (max_d_-min_d_)/(positions_.size()-1);
      drho  = (max_rho_-min_rho_)/(positions_.size()-1);
    }else{
      ddiam = 0.0;
      drho  = 0.0;
    }

    forAll(diameters_, i)
    {
        diameters_[i] = min_d_ + i*ddiam;
    }
    forAll(particleRhos_, i)
    {
        particleRhos_[i] = min_rho_ + i*drho;
    }

    // Determine volume of particles to inject
    this->volumeTotal_ = sum(pow3(diameters_))*pi/6.0;
}


template<class CloudType>
Foam::GGCInjection<CloudType>::GGCInjection
(
    const GGCInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    positionsFile_(im.positionsFile_),
    positions_(im.positions_),
    diameters_(im.diameters_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    U0_(im.U0_),
    min_d_(im.min_d_),
    max_d_(im.max_d_),
    min_rho_(im.min_rho_),
    max_rho_(im.max_rho_),
    particleRhos_(im.particleRhos_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::GGCInjection<CloudType>::~GGCInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::GGCInjection<CloudType>::timeEnd() const
{
    // Not used
    return this->SOI_;
}


template<class CloudType>
Foam::label Foam::GGCInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return positions_.size();
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::GGCInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    // All parcels introduced at SOI
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return this->volumeTotal_;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::GGCInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
)
{
    position = positions_[parcelI];
    cellOwner = injectorCells_[parcelI];
    tetFaceI = injectorTetFaces_[parcelI];
    tetPtI = injectorTetPts_[parcelI];
}


template<class CloudType>
void Foam::GGCInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = U0_;

    // set particle diameter
    parcel.d() = diameters_[parcelI];

    // set particle density
    parcel.rho() = particleRhos_[parcelI];
}


template<class CloudType>
bool Foam::GGCInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::GGCInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
