/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "PhaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
const Foam::wordList Foam::PhaseChangeModel<CloudType>::
enthalpyTransferTypeNames
(
    IStringStream
    (
        "("
            "latentHeat "
            "enthalpyDifference"
        ")"
    )()
);


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
typename Foam::PhaseChangeModel<CloudType>::enthalpyTransferType
Foam::PhaseChangeModel<CloudType>::wordToEnthalpyTransfer(const word& etName)
const
{
    forAll(enthalpyTransferTypeNames, i)
    {
        if (etName == enthalpyTransferTypeNames[i])
        {
            return enthalpyTransferType(i);
        }
    }

    FatalErrorIn
    (
        "PhaseChangeModel<CloudType>::enthalpyTransferType"
        "PhaseChangeModel<CloudType>::wordToEnthalpyTransfer(const word&) const"
    )   << "Unknown enthalpyType " << etName << ". Valid selections are:" << nl
        << enthalpyTransferTypeNames << exit(FatalError);

    return enthalpyTransferType(0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PhaseChangeModel<CloudType>::PhaseChangeModel
(
    CloudType& owner
)
:
    SubModelBase<CloudType>(owner),
    enthalpyTransfer_(etLatentHeat)
{}


template<class CloudType>
Foam::PhaseChangeModel<CloudType>::PhaseChangeModel
(
    const PhaseChangeModel<CloudType>& pcm
)
:
    SubModelBase<CloudType>(pcm),
    enthalpyTransfer_(pcm.enthalpyTransfer_)
{}


template<class CloudType>
Foam::PhaseChangeModel<CloudType>::PhaseChangeModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, type),
    enthalpyTransfer_
    (
        wordToEnthalpyTransfer(this->coeffDict().lookup("enthalpyTransfer"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PhaseChangeModel<CloudType>::~PhaseChangeModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const typename Foam::PhaseChangeModel<CloudType>::enthalpyTransferType&
Foam::PhaseChangeModel<CloudType>::enthalpyTransfer() const
{
    return enthalpyTransfer_;
}


template<class CloudType>
void Foam::PhaseChangeModel<CloudType>::calculate
(
    const scalar,
    const label,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    scalarField&
) const
{
    notImplemented
    (
        "void Foam::PhaseChangeModel<CloudType>::calculate"
        "("
            "const scalar, "
            "const label, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "const scalar, "
            "scalarField&"
        ") const"
    );
}


template<class CloudType>
Foam::scalar Foam::PhaseChangeModel<CloudType>::dh
(
    const label idc,
    const label idl,
    const label p,
    const label T
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PhaseChangeModelNew.C"

// ************************************************************************* //

