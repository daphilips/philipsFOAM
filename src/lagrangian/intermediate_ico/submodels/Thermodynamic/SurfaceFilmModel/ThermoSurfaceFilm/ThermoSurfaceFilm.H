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

Class
    Foam::ThermoSurfaceFilm

Description
    Thermo parcel surface film model.

    Responsible for:
    - injecting parcelss from the film model into the cloud, e.g. for dripping
    - parcel interaction with the film, e.g absorb, bounce, splash

    Splash model references:

        Bai and Gosman, `Mathematical modelling of wall films formed by
        impinging sprays', SAE 960626, 1996

        Bai et al, `Modelling off gasoline spray impingement', Atom. Sprays,
        vol 12, pp 1-27, 2002


SourceFiles
    ThermoSurfaceFilm.C
    ThermoSurfaceFilmI.H

\*---------------------------------------------------------------------------*/

#ifndef ThermoSurfaceFilm_H
#define ThermoSurfaceFilm_H

#include "SurfaceFilmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ThermoSurfaceFilm Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThermoSurfaceFilm
:
    public SurfaceFilmModel<CloudType>
{
public:

    // Public data

        //  Interaction type enumerations
        enum interactionType
        {
            itAbsorb,
            itBounce,
            itSplashBai
        };

        //- Word descriptions of interaction type names
        static wordList interactionTypeNames_;


    // Public Member Functions

        // Return interaction type enum from word
        interactionType interactionTypeEnum(const word& it) const;

        // Return word from interaction type enum
        word interactionTypeStr(const interactionType& it) const;


protected:

    // Protected data

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

        //- Reference to the cloud random number generator
        cachedRandom& rndGen_;

        //- Reference to the cloud thermo package
        const SLGThermo& thermo_;


        // Cached injector fields per film patch

            //- Film temperature / patch face
            scalarList TFilmPatch_;

            //- Film specific heat capacity / patch face
            scalarList CpFilmPatch_;


        // Interaction model data

            //- Interaction type enumeration
            interactionType interactionType_;

            //- Film thickness beyond which patch is assumed to be wet
            scalar deltaWet_;

            //- Splash parcel type label - id assigned to identify parcel for
            //  post-processing. If not specified, defaults to originating cloud
            //  type
            label splashParcelType_;

            //- Number of new parcels resulting from splash event
            label parcelsPerSplash_;


            // Surface roughness coefficient typically in the range 1300 - 5200
            // and decreases with increasing surface roughness

                //- Dry surface roughness coefficient
                //  = 2630 for dry interaction (ref. Bai)
                scalar Adry_;

                //- Wet surface roughness coefficient
                //  = 1320 for wet interaction (ref. Bai)
                scalar Awet_;


            //- Skin friction typically in the range 0.6 < Cf < 0.8
            scalar Cf_;

            //- Counter for number of new splash parcels
            label nParcelsSplashed_;


    // Protected Member Functions

        //- Return a vector tangential to input vector, v
        vector tangentVector(const vector& v) const;

        //- Return splashed parcel direction
        vector splashDirection
        (
            const vector& tanVec1,
            const vector& tanVec2,
            const vector& nf
        ) const;


        // Interaction models

            //- Absorb parcel into film
            void absorbInteraction
            (
                regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
                const parcelType& p,
                const polyPatch& pp,
                const label faceI,
                const scalar mass,
                bool& keepParticle
            );

            //- Bounce parcel (flip parcel normal velocity)
            void bounceInteraction
            (
                parcelType& p,
                const polyPatch& pp,
                const label faceI,
                bool& keepParticle
            ) const;

            //- Parcel interaction with dry surface
            void drySplashInteraction
            (
                regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
                const parcelType& p,
                const polyPatch& pp,
                const label faceI,
                bool& keepParticle
            );

            //- Parcel interaction with wetted surface
            void wetSplashInteraction
            (
                regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
                parcelType& p,
                const polyPatch& pp,
                const label faceI,
                bool& keepParticle
            );

            //- Bai parcel splash interaction model
            void splashInteraction
            (
                regionModels::surfaceFilmModels::surfaceFilmModel& filmModel,
                const parcelType& p,
                const polyPatch& pp,
                const label faceI,
                const scalar mRatio,
                const scalar We,
                const scalar Wec,
                const scalar sigma,
                bool& keepParticle
            );


        // Injection from sheet (ejection) helper functions

            //- Cache the film fields in preparation for injection
            virtual void cacheFilmFields
            (
                const label filmPatchI,
                const label primaryPatchI,
                const mapDistribute& distMap,
                const regionModels::surfaceFilmModels::surfaceFilmModel&
                    filmModel
            );

            //- Set the individual parcel properties
            virtual void setParcelProperties
            (
                parcelType& p,
                const label filmFaceI
            ) const;


public:

    //- Runtime type information
    TypeName("thermoSurfaceFilm");


    // Constructors

        //- Construct from components
        ThermoSurfaceFilm
        (
            const dictionary& dict,
            CloudType& owner,
            const dimensionedVector& g
        );

        //- Construct copy
        ThermoSurfaceFilm(const ThermoSurfaceFilm<CloudType>& sfm);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<SurfaceFilmModel<CloudType> > clone()
        {
            return autoPtr<SurfaceFilmModel<CloudType> >
            (
                new ThermoSurfaceFilm<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ThermoSurfaceFilm();


    // Member Functions

        // Evaluation

            //- Transfer parcel from cloud to surface film
            //  Returns true if parcel is to be transferred
            virtual bool transferParcel
            (
                parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );


        // I-O

            //- Write surface film info to stream
            virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ThermoSurfaceFilm.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
