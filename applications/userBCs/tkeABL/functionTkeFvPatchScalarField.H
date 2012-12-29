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
    Boundary condition specifies a function Turbulence kinetic energy
    (profile for ABL modelling), given roughness height, constants 
    A and B, and the height coordinate y

SourceFiles
    functionVelocityFvPatchVectorField.C


\*---------------------------------------------------------------------------*/

#ifndef functionTkeFvPatchScalarField_H
#define functionTkeFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class functionTkeFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class functionTkeFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data
        
        //- roughness height y0
        scalar yZero_;
        
        //- constant A;
        scalar A_;

        //- constant B;
        scalar B_;
        
        //- Direction of the y-coordinate
        vector y_;


public:

    //- Runtime type information
    TypeName("functionTke");


    // Constructors

        //- Construct from patch and internal field
        functionTkeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        functionTkeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given functionTkeFvPatchScalarField
        //  onto a new patch
        functionTkeFvPatchScalarField
        (
            const functionTkeFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new functionTkeFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        functionTkeFvPatchScalarField
        (
            const functionTkeFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new functionTkeFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        
        //- Return roughness height
        scalar& yZero()
        {
            return yZero_;
        }

        //- Return constant A
        scalar& A()
        {
            return A_;
        }

        //- Return constant B
        scalar& B()
        {
            return B_;
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
