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

Class
    Foam::StandardWallInteraction

Description
    Wall interaction model. Three choices:
    - rebound - optionally specify elasticity and resitution coefficients
    - stick   - particles assigined zero velocity
    - escape  - remove particle from the domain

    Example usage:

        StandardWallInteractionCoeffs
        {
            type        rebound; // stick, escape
            e           1;       // optional - elasticity coeff
            mu          0;       // optional - restitution coeff
        }

\*---------------------------------------------------------------------------*/

#ifndef StandardWallInteraction_H
#define StandardWallInteraction_H

#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                   Class StandardWallInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class StandardWallInteraction
:
    public PatchInteractionModel<CloudType>
{
protected:

    // Protected data

        //- Interaction type
        typename PatchInteractionModel<CloudType>::interactionType
            interactionType_;

        //- Elasticity coefficient
        scalar e_;

        //- Restitution coefficient
        scalar mu_;


        // Counters for initial particle fates

            //- Number of parcels escaped
            label nEscape0_;

            //- Mass of parcels escaped
            scalar massEscape0_;

            //- Number of parcels stuck to patches
            label nStick0_;

            //- Mass of parcels stuck to patches
            scalar massStick0_;


        // Counters for particle fates

            //- Number of parcels escaped
            label nEscape_;

            //- Mass of parcels escaped
            scalar massEscape_;

            //- Number of parcels stuck to patches
            label nStick_;

            //- Mass of parcels stuck to patches
            scalar massStick_;


    // Protected Member Functions

        //- Read interaction properties from file
        void readProps();

        //- Write interaction properties to file
        void writeProps
        (
            const label nEscape,
            const scalar massEscape,
            const label nStick,
            const scalar massStick
        ) const;


public:

    //- Runtime type information
    TypeName("standardWallInteraction");


    // Constructors

        //- Construct from dictionary
        StandardWallInteraction(const dictionary& dict, CloudType& cloud);

        //- Construct copy from owner cloud and patch interaction model
        StandardWallInteraction(const StandardWallInteraction<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType> > clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType> >
            (
                new StandardWallInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~StandardWallInteraction();


    // Member Functions

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle,
            const scalar trackFraction,
            const tetIndices& tetIs
        );


        // I-O

            //- Write patch interaction info to stream
            virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "StandardWallInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
