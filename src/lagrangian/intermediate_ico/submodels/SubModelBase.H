/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2011 OpenCFD Ltd.
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
    Foam::SubModelBase

Description
    Base class for cloud sub-models

SourceFiles
    SubModelBase.C

\*---------------------------------------------------------------------------*/

#ifndef SubModelBase_H
#define SubModelBase_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CloudType>
class SubModelBase;

/*---------------------------------------------------------------------------*\
                       Class SubModelBase Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SubModelBase
{
protected:

    // Protected Data

        //- Reference to the cloud
        CloudType& owner_;

        //- Reference to the cloud dictionary
        const dictionary& dict_;

        //- Coefficients dictionary
        const dictionary& coeffDict_;


public:

    // Constructors

        //- Construct null from owner cloud
        SubModelBase(CloudType& owner);

        //- Construct from owner cloud, dictionary, and model type name
        SubModelBase
        (
            CloudType& owner,
            const dictionary& dict,
            const word& name,
            const word& dictExt = "Coeffs"
        );

        //- Construct as copy
        SubModelBase(const SubModelBase<CloudType>& smb);


    //- Destructor
    virtual ~SubModelBase();

    //- Type of cloud this model was instantiated for
    typedef CloudType cloudType;


    // Member Functions

        // Access

            //- Return const access to the cloud dictionary
            const dictionary& dict() const;

            //- Return const access to the owner cloud
            const CloudType& owner() const;

            //- Return const access to the coefficients dictionary
            const dictionary& coeffDict() const;

            //- Returns true if defaultCoeffs is true and outputs on printMsg
            bool defaultCoeffs(const bool printMsg) const;

            //- Return the model 'active' status - default active = true
            virtual bool active() const;

            //- Cache dependant sub-model fields
            virtual void cacheFields(const bool store);


        // Edit

            //- Return non-const access to the owner cloud for manipulation
            CloudType& owner();


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SubModelBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
