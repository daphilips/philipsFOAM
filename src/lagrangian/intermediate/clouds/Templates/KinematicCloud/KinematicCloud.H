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
    Foam::KinematicCloud

Description
    Templated base class for kinematic cloud

    - cloud function objects

    - particle forces, e.g.
      - buoyancy
      - drag
      - pressure gradient
      - ...

    - sub-models:
      - dispersion model
      - injection model
      - patch interaction model
      - surface film model

SourceFiles
    KinematicCloudI.H
    KinematicCloud.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicCloud_H
#define KinematicCloud_H

#include "particle.H"
#include "Cloud.H"
#include "kinematicCloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "cachedRandom.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "IntegrationSchemesFwd.H"
#include "cloudSolution.H"

#include "ParticleForceList.H"
#include "CloudFunctionObjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class DispersionModel;

template<class CloudType>
class InjectionModel;

template<class CloudType>
class PatchInteractionModel;

template<class CloudType>
class SurfaceFilmModel;


/*---------------------------------------------------------------------------*\
                       Class KinematicCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class KinematicCloud
:
    public CloudType,
    public kinematicCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef KinematicCloud<CloudType> kinematicCloudType;

        //- Force models type
        typedef ParticleForceList<KinematicCloud<CloudType> > forceType;

        //- Function object type
        typedef CloudFunctionObjectList<KinematicCloud<CloudType> >
            functionType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<KinematicCloud<CloudType> > cloudCopyPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        KinematicCloud(const KinematicCloud&);

        //- Disallow default bitwise assignment
        void operator=(const KinematicCloud&);


protected:

    // Protected data

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- Solution properties
        cloudSolution solution_;

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;

        //- Sub-models dictionary
        const dictionary subModelProperties_;

        //- Random number generator - used by some injection routines
        cachedRandom rndGen_;

        //- Cell occupancy information for each parcel, (demand driven)
        autoPtr<List<DynamicList<parcelType*> > > cellOccupancyPtr_;


        // References to the carrier gas fields

            //- Density [kg/m3]
            const volScalarField& rho_;

            //- Velocity [m/s]
            const volVectorField& U_;

            //- Dynamic viscosity [Pa.s]
            const volScalarField& mu_;


        // Environmental properties

            //- Gravity
            const dimensionedVector& g_;

            //- Averaged ambient domain pressure
            scalar pAmbient_;


        //- Optional particle forces
        forceType forces_;

        //- Optional cloud function objects
        functionType functions_;


        // References to the cloud sub-models

            //- Dispersion model
            autoPtr<DispersionModel<KinematicCloud<CloudType> > >
                dispersionModel_;

            //- Injector model
            autoPtr<InjectionModel<KinematicCloud<CloudType> > >
                injectionModel_;

            //- Patch interaction model
            autoPtr<PatchInteractionModel<KinematicCloud<CloudType> > >
                patchInteractionModel_;

            //- Surface film model
            autoPtr<SurfaceFilmModel<KinematicCloud<CloudType> > >
                surfaceFilmModel_;


        // Reference to the particle integration schemes

            //- Velocity integration
            autoPtr<vectorIntegrationScheme> UIntegrator_;


        // Sources

            //- Momentum
            autoPtr<DimensionedField<vector, volMesh> > UTrans_;

            //- Coefficient for carrier phase U equation
            autoPtr<DimensionedField<scalar, volMesh> > UCoeff_;


        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Solve the cloud - calls all evolution functions
            template<class TrackData>
            void solve(TrackData& td);

            //- Build the cellOccupancy
            void buildCellOccupancy();

            //- Update (i.e. build) the cellOccupancy if it has
            //  already been used
            void updateCellOccupancy();

            //- Evolve the cloud
            template<class TrackData>
            void evolveCloud(TrackData& td);

            //- Post-evolve
            void postEvolve();

            //- Reset state of cloud
            void cloudReset(KinematicCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        KinematicCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        KinematicCloud
        (
            KinematicCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        KinematicCloud
        (
            const fvMesh& mesh,
            const word& name,
            const KinematicCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<parcelType> >
            (
                new KinematicCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType> >
            (
                new KinematicCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~KinematicCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const KinematicCloud& cloudCopy() const;

            //- Switch to specify if particles of the cloud can return
            //  non-zero wall distance values - true for kinematic parcels
            virtual bool hasWallImpactDistance() const;


            // References to the mesh and databases

                //- Return reference to the mesh
                inline const fvMesh& mesh() const;

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

                //- Return const access to the solution properties
                inline const cloudSolution& solution() const;

                //- Return access to the solution properties
                inline cloudSolution& solution();

                //- Return the constant properties
                inline const typename parcelType::constantProperties&
                    constProps() const;

                    //- Return reference to the sub-models dictionary
                inline const dictionary& subModelProperties() const;


            // Cloud data

                //- Return reference to the random object
                inline cachedRandom& rndGen();

                //- Return the cell occupancy information for each
                //  parcel, non-const access, the caller is
                //  responsible for updating it for its own purposes
                //  if particles are removed or created.
                inline List<DynamicList<parcelType*> >& cellOccupancy();


            // References to the carrier gas fields

                //- Return carrier gas velocity
                inline const volVectorField& U() const;

                //- Return carrier gas density
                inline const volScalarField& rho() const;

                //- Return carrier gas dynamic viscosity
                inline const volScalarField& mu() const;


            // Environmental properties

                //- Gravity
                inline const dimensionedVector& g() const;

                //- Return const-access to the ambient pressure
                inline scalar pAmbient() const;

                //- Return reference to the ambient pressure
                inline scalar& pAmbient();


            //- Optional particle forces
//            inline const typename parcelType::forceType& forces() const;
            inline const forceType& forces() const;

            //- Optional cloud function objects
            inline functionType& functions();


            // Sub-models

                //- Return const-access to the dispersion model
                inline const DispersionModel<KinematicCloud<CloudType> >&
                    dispersion() const;

                //- Return reference to the dispersion model
                inline DispersionModel<KinematicCloud<CloudType> >&
                    dispersion();

                //- Return const access to the injection model
                inline const InjectionModel<KinematicCloud<CloudType> >&
                    injection() const;

                //- Return reference to the injection model
                inline InjectionModel<KinematicCloud<CloudType> >&
                    injection();

                //- Return const-access to the patch interaction model
                inline const PatchInteractionModel<KinematicCloud<CloudType> >&
                    patchInteraction() const;

                //- Return reference to the patch interaction model
                inline PatchInteractionModel<KinematicCloud<CloudType> >&
                    patchInteraction();

                //- Return const-access to the surface film model
                inline const SurfaceFilmModel<KinematicCloud<CloudType> >&
                    surfaceFilm() const;

                //- Return reference to the surface film model
                inline SurfaceFilmModel<KinematicCloud<CloudType> >&
                    surfaceFilm();


            // Integration schemes

                //-Return reference to velocity integration
                inline const vectorIntegrationScheme& UIntegrator() const;


            // Sources

                // Momentum

                    //- Return reference to momentum source
                    inline DimensionedField<vector, volMesh>& UTrans();

                    //- Return const reference to momentum source
                    inline const DimensionedField<vector, volMesh>&
                        UTrans() const;

                     //- Return coefficient for carrier phase U equation
                    inline DimensionedField<scalar, volMesh>& UCoeff();

                    //- Return const coefficient for carrier phase U equation
                    inline const DimensionedField<scalar, volMesh>&
                        UCoeff() const;

                    //- Return tmp momentum source term
                    inline tmp<fvVectorMatrix> SU(volVectorField& U) const;


        // Check

            //- Total mass injected
            inline scalar massInjected() const;

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            //- Total rotational kinetic energy in the system
            inline scalar rotationalKineticEnergyOfSystem() const;

            //- Penetration for percentage of the current total mass
            inline scalar penetration(const scalar& prc) const;

            //- Mean diameter Dij
            inline scalar Dij(const label i, const label j) const;


            // Fields

                //- Return the particle volume fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> theta() const;

                //- Return the particle mass fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> alpha() const;

                //- Return the particle effective density field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> rhoEff() const;


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Relax field
            template<class Type>
            void relax
            (
                DimensionedField<Type, volMesh>& field,
                const DimensionedField<Type, volMesh>& field0,
                const word& name
            ) const;

            //- Scale field
            template<class Type>
            void scale
            (
                DimensionedField<Type, volMesh>& field,
                const word& name
            ) const;

            //- Apply relaxation to (steady state) cloud sources
            void relaxSources(const KinematicCloud<CloudType>& cloudOldTime);

            //- Apply scaling to (transient) cloud sources
            void scaleSources();

            //- Pre-evolve
            void preEvolve();

            //- Evolve the cloud
            void evolve();

            //- Particle motion
            template<class TrackData>
            void motion(TrackData& td);

            //- Print cloud information
            void info() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KinematicCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "KinematicCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
