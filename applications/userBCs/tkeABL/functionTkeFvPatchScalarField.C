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

\*---------------------------------------------------------------------------*/

#include "functionTkeFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

functionTkeFvPatchScalarField::functionTkeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    yZero_(0.00075),
    A_(0.025),
    B_(0.41),
    y_(0, 0, 1)
{}


functionTkeFvPatchScalarField::functionTkeFvPatchScalarField
(
    const functionTkeFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    yZero_(ptf.yZero_),
    A_(ptf.A_),
    B_(ptf.B_),
    y_(ptf.y_)
{}


functionTkeFvPatchScalarField::functionTkeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    yZero_(readScalar(dict.lookup("yZero"))),
    A_(readScalar(dict.lookup("A"))),
    B_(readScalar(dict.lookup("B"))),
    y_(dict.lookup("y"))
{
    if (mag(y_) < SMALL)
    {
        FatalErrorIn("functionTkeFvPatchScalarField(dict)")
            << "y given with zero size not correct"
            << abort(FatalError);
    }

    y_ /= mag(y_);

    evaluate();
}


functionTkeFvPatchScalarField::functionTkeFvPatchScalarField
(
    const functionTkeFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    yZero_(fcvpvf.yZero_),
    A_(fcvpvf.A_),
    B_(fcvpvf.B_),
    y_(fcvpvf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void functionTkeFvPatchScalarField::updateCoeffs()
{

//    Info << "Function called." << endl;

    // Get range and orientation
    //boundBox bb(patch().patch().localPoints(), false);
    //vector ctr = 0.5*(bb.max() + bb.min());

    const vectorField& c = patch().Cf();

    // Calculate argument of square root
    scalarField coord = (A_*log((c & y_)+yZero_)+B_);
        
    //scalarField::operator=sqrt(coord);
    scalarField::operator = (sqrt(coord));
}


// Write
void functionTkeFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("yZero")
        << yZero_ << token::END_STATEMENT << nl;
    os.writeKeyword("A")
        << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("B")
        << B_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, functionTkeFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
