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
    Foam::ReactingParcel

Description
    Reacting parcel class with one/two-way coupling with the continuous
    phase.

SourceFiles
    ReactingParcelI.H
    ReactingParcel.C
    ReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingParcel_H
#define ReactingParcel_H

#include "particle.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                        Class ReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingParcel
:
    public ParcelType
{
public:

    //- Class to hold reacting particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private data

            //- Minimum pressure [Pa]
            scalar pMin_;

            //- Constant volume flag - e.g. during mass transfer
            Switch constantVolume_;

            //- Vaporisation temperature [K]
            scalar Tvap_;

            //- Boiling point [K]
            scalar Tbp_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Constructor from dictionary
            constantProperties
            (
                const dictionary& parentDict,
                const bool readFields = true
            );

            //- Construct from components
            constantProperties
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
            );


        // Access

            //- Return const access to the minimum pressure
            inline scalar pMin() const;

            //- Return const access to the constant volume flag
            inline Switch constantVolume() const;

            //- Return const access to the vaporisation temperature
            inline scalar Tvap() const;

            //- Return const access to the boiling point
            inline scalar Tbp() const;
    };


    template<class CloudType>
    class TrackingData
    :
        public ParcelType::template TrackingData<CloudType>
    {
    private:

        // Private data

            // Interpolators for continuous phase fields

                //- Interpolator for continuous phase pressure field
                autoPtr<interpolation<scalar> > pInterp_;


    public:

        typedef typename ParcelType::template TrackingData<CloudType>::trackPart
            trackPart;

        // Constructors

            //- Construct from components
            inline TrackingData
            (
                CloudType& cloud,
                trackPart part = ParcelType::template
                    TrackingData<CloudType>::tpLinearTrack
            );


        // Member functions

            //- Return const access to the interpolator for continuous phase
            //  pressure field
            inline const interpolation<scalar>& pInterp() const;
    };


protected:

    // Protected data

        // Parcel properties

            //- Initial particle mass [kg]
            scalar mass0_;

            //- Mass fractions of mixture []
            scalarField Y_;


        // Cell-based quantities

            //- Pressure [Pa]
            scalar pc_;


    // Protected Member Functions

        //- Calculate Phase change
        template<class TrackData>
        void calcPhaseChange
        (
            TrackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar Re,           // Reynolds number
            const scalar Ts,           // Surface temperature
            const scalar nus,          // Surface kinematic viscosity
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label idPhase,       // id of phase involved in phase change
            const scalar YPhase,       // total mass fraction
            const scalarField& YComponents, // component mass fractions
            scalarField& dMassPC,      // mass transfer - local to particle
            scalar& Sh,                // explicit particle enthalpy source
            scalar& N,                 // flux of species emitted from particle
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        );

        //- Update mass fraction
        scalar updateMassFraction
        (
            const scalar mass0,
            const scalarField& dMass,
            scalarField& Y
        ) const;


public:

    // Static data members

        //- String representation of properties
        static string propHeader;

        //- Runtime type information
        TypeName("ReactingParcel");


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline ReactingParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components
        inline ReactingParcel
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
        );

        //- Construct from Istream
        ReactingParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ReactingParcel
        (
            const ReactingParcel& p,
            const polyMesh& mesh
        );

        //- Construct as a copy
        ReactingParcel(const ReactingParcel& p);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ReactingParcel<ParcelType>(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>
            (
                new ReactingParcel<ParcelType>(*this, mesh)
            );
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ReactingParcel<ParcelType> > operator()(Istream& is) const
            {
                return autoPtr<ReactingParcel<ParcelType> >
                (
                    new ReactingParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to initial mass
            inline scalar mass0() const;

            //- Return const access to mass fractions of mixture
            inline const scalarField& Y() const;

            //- Return the owner cell pressure
            inline scalar pc() const;

            //- Return reference to the owner cell pressure
            inline scalar& pc();


        // Edit

            //- Return access to initial mass
            inline scalar& mass0();

            //- Return access to mass fractions of mixture
            inline scalarField& Y();


        // Main calculation loop

            //- Set cell values
            template<class TrackData>
            void setCellValues
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct cell values using latest transfer information
            template<class TrackData>
            void cellValueSourceCorrection
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct surface values due to emitted species
            template<class TrackData>
            void correctSurfaceValues
            (
                TrackData& td,
                const label cellI,
                const scalar T,
                const scalarField& Cs,
                scalar& rhos,
                scalar& mus,
                scalar& Prs,
                scalar& kappas
            );

            //- Update parcel properties over the time interval
            template<class TrackData>
            void calc
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );


        // I-O

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Write - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingParcelI.H"
#include "ReactingParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
