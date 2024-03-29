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
    Foam::SurfaceReactionModel

Description
    Templated surface reaction model class

SourceFiles
    SurfaceReactionModel.C
    SurfaceReactionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef SurfaceReactionModel_H
#define SurfaceReactionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "SubModelBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class SurfaceReactionModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SurfaceReactionModel
:
    public SubModelBase<CloudType>
{
public:

    //-Runtime type information
    TypeName("surfaceReactionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        SurfaceReactionModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        SurfaceReactionModel(CloudType& owner);

        //- Construct from dictionary
        SurfaceReactionModel
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& type
        );

        //- Construct copy
        SurfaceReactionModel(const SurfaceReactionModel<CloudType>& srm);

        //- Construct and return a clone
        virtual autoPtr<SurfaceReactionModel<CloudType> > clone() const
        {
            return autoPtr<SurfaceReactionModel<CloudType> >
            (
                new SurfaceReactionModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~SurfaceReactionModel();


    //- Selector
    static autoPtr<SurfaceReactionModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& cloud
    );


    // Member Functions

        //- Update surface reactions
        //  Returns the heat of reaction
        virtual scalar calculate
        (
            const scalar dt,
            const label cellI,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalar N,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSurfaceReactionModel(CloudType)                                   \
                                                                              \
    typedef CloudType::reactingMultiphaseCloudType                            \
        reactingMultiphaseCloudType;                                          \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SurfaceReactionModel<reactingMultiphaseCloudType>,                    \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        SurfaceReactionModel<reactingMultiphaseCloudType>,                    \
        dictionary                                                            \
    );


#define makeSurfaceReactionModelType(SS, CloudType)                           \
                                                                              \
    typedef CloudType::reactingMultiphaseCloudType                            \
        reactingMultiphaseCloudType;                                          \
    defineNamedTemplateTypeNameAndDebug(SS<reactingMultiphaseCloudType>, 0);  \
                                                                              \
    SurfaceReactionModel<reactingMultiphaseCloudType>::                       \
        adddictionaryConstructorToTable<SS<reactingMultiphaseCloudType> >     \
        add##SS##CloudType##reactingMultiphaseCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SurfaceReactionModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
