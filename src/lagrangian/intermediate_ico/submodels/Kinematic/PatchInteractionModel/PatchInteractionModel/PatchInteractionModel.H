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
    Foam::PatchInteractionModel

Description
    Templated patch interaction model class

SourceFiles
    PatchInteractionModel.C
    PatchInteractionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PatchInteractionModel_H
#define PatchInteractionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "polyPatch.H"
#include "wallPolyPatch.H"
#include "tetIndices.H"
#include "SubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class PatchInteractionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchInteractionModel
:
    public SubModelBase<CloudType>
{
public:

    // Public enumerations

        // Interaction types
        enum interactionType
        {
            itRebound,
            itStick,
            itEscape,
            itOther
        };

        static wordList interactionTypeNames_;


private:

    // Private data

        //- Name of velocity field - default = "U"
        const word UName_;


public:

    //- Runtime type information
    TypeName("patchInteractionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PatchInteractionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        PatchInteractionModel(CloudType& owner);

        //- Construct from components
        PatchInteractionModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        PatchInteractionModel(const PatchInteractionModel<CloudType>& pim);

        //- Construct and return a clone
        virtual autoPtr<PatchInteractionModel<CloudType> > clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType> >
            (
                new PatchInteractionModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PatchInteractionModel();


    //- Selector
    static autoPtr<PatchInteractionModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        //- Return name of velocity field
        const word& UName() const;


    // Member Functions

        //- Convert interaction result to word
        static word interactionTypeToWord(const interactionType& itEnum);

        //- Convert word to interaction result
        static interactionType wordToInteractionType(const word& itWord);

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

        //- Calculate the patch normal and velocity to interact with,
        //  accounting for patch motion if required.
        void patchData
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            const scalar trackFraction,
            const tetIndices& tetIs,
            vector& normal,
            vector& Up
        ) const;


        // I-O

            //- Write patch interaction info to stream
            virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePatchInteractionModel(CloudType)                                  \
                                                                              \
    typedef CloudType::kinematicCloudType kinematicCloudType;                 \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        PatchInteractionModel<kinematicCloudType>,                            \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        PatchInteractionModel<kinematicCloudType>,                            \
        dictionary                                                            \
    );


#define makePatchInteractionModelType(SS, CloudType)                          \
                                                                              \
    typedef CloudType::kinematicCloudType kinematicCloudType;                 \
    defineNamedTemplateTypeNameAndDebug(SS<kinematicCloudType>, 0);           \
                                                                              \
    PatchInteractionModel<kinematicCloudType>::                               \
        adddictionaryConstructorToTable<SS<kinematicCloudType> >              \
            add##SS##CloudType##kinematicCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PatchInteractionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
