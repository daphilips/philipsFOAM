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

#ifndef functionOmegaFvPatchScalarField_H
#define functionOmegaFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class functionOmegaFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class functionOmegaFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- friction velocity uStar
        scalar uStar_;

        //- von Karman constant
        scalar kappa_;

        //- Cmu
        scalar Cmu_;
        
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
    TypeName("functionOmega");


    // Constructors

        //- Construct from patch and internal field
        functionOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        functionOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given functionOmegaFvPatchScalarField
        //  onto a new patch
        functionOmegaFvPatchScalarField
        (
            const functionOmegaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new functionOmegaFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        functionOmegaFvPatchScalarField
        (
            const functionOmegaFvPatchScalarField&,
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
                new functionOmegaFvPatchScalarField(*this, iF)
            );
        }


    // Member functions


        //- Return roughness height
        scalar& uStar()
        {
            return uStar_;
        }

        //- Return roughness height
        scalar& Cmu()
        {
            return Cmu_;
        }

        //- Return roughness height
        scalar& kappa()
        {
            return kappa_;
        }
        
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
