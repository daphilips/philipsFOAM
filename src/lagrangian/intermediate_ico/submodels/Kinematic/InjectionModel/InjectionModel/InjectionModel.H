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
    Foam::InjectionModel

Description
    Templated injection model class.

    The injection model nominally describes the parcel:
    - position
    - diameter
    - velocity
    In this case, the fullyDescribed() flag should be set to 0 (false). When
    the parcel is then added to the cloud, the remaining properties are
    populated using values supplied in the constant properties.

    If, however, all of a parcel's properties are described in the model, the
    fullDescribed() flag should be set to 1 (true).


SourceFiles
    InjectionModel.C
    InjectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef InjectionModel_H
#define InjectionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "SubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class InjectionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InjectionModel
:
    public SubModelBase<CloudType>
{
public:

    //- Convenience typedef for parcelType
    typedef typename CloudType::parcelType parcelType;


    // Enumerations

        //- Parcel basis representation options
        //  i.e constant number of particles OR constant mass per parcel
        enum parcelBasis
        {
            pbNumber,
            pbMass,
            pbFixed
        };


private:

    // Private Member Functions

        //- Read injector properties from previous run (if applicable)
        void readProps();

        //- Write injector properties
        void writeProps();


protected:

    // Protected data

        // Global injection properties

            //- Start of injection [s]
            scalar SOI_;

            //- Total volume of particles introduced by this injector [m^3]
            //  - scaled to ensure massTotal is achieved
            scalar volumeTotal_;

            //- Total mass to inject [kg]
            scalar massTotal_;

            //- Total mass injected to date [kg]
            scalar massInjected_;


        // Counters

            //- Number of injections counter
            label nInjections_;

            //- Running counter of total number of parcels added
            label parcelsAddedTotal_;


        // Injection properties per Lagrangian time step

            //- Parcel basis enumeration
            parcelBasis parcelBasis_;

            //- nParticle to assign to parcels when the 'fixed' basis
            //  is selected
            scalar nParticleFixed_;

            //- Continuous phase time at start of injection time step [s]
            scalar time0_;

            //- Time at start of injection time step [s]
            scalar timeStep0_;


    // Protected Member Functions

        //- Additional flag to identify whether or not injection of parcelI is
        //  permitted
        virtual bool validInjection(const label parcelI);

        //- Determine properties for next time step/injection interval
        virtual void prepareForNextTimeStep
        (
            const scalar time,
            label& newParcels,
            scalar& newVolume
        );

        //- Find the cell that contains the supplied position
        //  Will modify position slightly towards the owner cell centroid to
        //  ensure that it lies in a cell and not edge/face
        virtual bool findCellAtPosition
        (
            label& cellI,
            label& tetFaceI,
            label& tetPtI,
            vector& position,
            bool errorOnNotFound = true
        );

        //- Set number of particles to inject given parcel properties
        virtual scalar setNumberOfParticles
        (
            const label parcels,
            const scalar volume,
            const scalar diameter,
            const scalar rho
        );

        //- Post injection checks
        virtual void postInjectCheck
        (
            const label parcelsAdded,
            const scalar massAdded
        );


public:

    //- Runtime type information
    TypeName("injectionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        InjectionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        InjectionModel(CloudType& owner);

        //- Construct from dictionary
        InjectionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        InjectionModel(const InjectionModel<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType> > clone() const
        {
            return autoPtr<InjectionModel<CloudType> >
            (
                new InjectionModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~InjectionModel();


    //- Selector
    static autoPtr<InjectionModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // Global information

            //- Return the start-of-injection time
            inline scalar timeStart() const;

            //- Return the total volume to be injected across the event
            inline scalar volumeTotal() const;

            //- Return mass of particles to introduce
            inline scalar massTotal() const;

            //- Return mass of particles injected (cumulative)
            inline scalar massInjected() const;

            //- Return the end-of-injection time
            virtual scalar timeEnd() const;

            //- Number of parcels to introduce relative to SOI
            virtual label parcelsToInject
            (
                const scalar time0,
                const scalar time1
            );

            //- Volume of parcels to introduce relative to SOI
            virtual scalar volumeToInject
            (
                const scalar time0,
                const scalar time1
            );

            //- Return the average parcel mass over the injection period
            scalar averageParcelMass();


            // Counters

                //- Return the number of injections
                inline label nInjections() const;

                //- Return the total number parcels added
                inline label parcelsAddedTotal() const;


        // Per-injection event functions

            //- Main injection loop
            template<class TrackData>
            void inject(TrackData& td);

            //- Main injection loop - steady-state
            template<class TrackData>
            void injectSteadyState(TrackData& td, const scalar trackTime);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFaceI,
                label& tetPtI
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInjectionModel(CloudType)                                         \
                                                                              \
    typedef CloudType::kinematicCloudType kinematicCloudType;                 \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        InjectionModel<kinematicCloudType>,                                   \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        InjectionModel<kinematicCloudType>,                                   \
        dictionary                                                            \
    );


#define makeInjectionModelType(SS, CloudType)                                 \
                                                                              \
    typedef CloudType::kinematicCloudType kinematicCloudType;                 \
    defineNamedTemplateTypeNameAndDebug(SS<kinematicCloudType>, 0);           \
                                                                              \
    InjectionModel<kinematicCloudType>::                                      \
        adddictionaryConstructorToTable<SS<kinematicCloudType> >              \
        add##SS##CloudType##kinematicCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InjectionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InjectionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
