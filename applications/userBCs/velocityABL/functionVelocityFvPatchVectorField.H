/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    functionVelocityFvPatchVectorField

Description
    Boundary condition specifies a function velocity inlet profile
    (logarithmic profile for ABL modelling), given friction velocity,
    von Karman constant, roughness height, flow direction n and 
    direction of the height coordinate y

SourceFiles
    functionVelocityFvPatchVectorField.C


\*---------------------------------------------------------------------------*/

#ifndef functionVelocityFvPatchVectorField_H
#define functionVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class functionVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class functionVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Friction Velocity
        scalar uStar_;
        
        //- von Karman constant
        scalar kappa_;
        
        //- roughness height y0
        scalar yZero_;
        
        //- Flow direction
        vector n_;

        //- Direction of the y-coordinate
        vector y_;


public:

    //- Runtime type information
    TypeName("logarithmicVelocity");


    // Constructors

        //- Construct from patch and internal field
        functionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        functionVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given functionVelocityFvPatchVectorField
        //  onto a new patch
        functionVelocityFvPatchVectorField
        (
            const functionVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new functionVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        functionVelocityFvPatchVectorField
        (
            const functionVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new functionVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        
        //- Return friction velocity
        scalar& uStar()
        {
            return uStar_;
        }
        
        //- Return von karman constant
        scalar& kappa()
        {
            return kappa_;
        }
        
        //- Return roughness height
        scalar& yZero()
        {
            return yZero_;
        }

        //- Return flow direction
        vector& n()
        {
            return n_;
        }

        //- Return y direction
        vector& y()
        {
            return y_;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
