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

#include "functionVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

functionVelocityFvPatchVectorField::functionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(0),
    n_(1, 0, 0),
    y_(0, 0, 1)
{}


functionVelocityFvPatchVectorField::functionVelocityFvPatchVectorField
(
    const functionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_)
{}


functionVelocityFvPatchVectorField::functionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(readScalar(dict.lookup("maxValue"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("functionVelocityFvPatchVectorField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);

    evaluate();
}


functionVelocityFvPatchVectorField::functionVelocityFvPatchVectorField
(
    const functionVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    maxValue_(fcvpvf.maxValue_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void functionVelocityFvPatchVectorField::updateCoeffs()
{

    Info << "Function called." << endl;

    // Get range and orientation
    //boundBox bb(patch().patch().localPoints(), false);

    //vector ctr = 0.5*(bb.max() + bb.min());

    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate for the parabolic profile
    //scalarField coord = 2*((c - ctr) & y_)/((bb.max() - bb.min()) & y_);
    scalarField coord = ((0.00075 + (c & y_))/0.00075);
    
    //vectorField::operator=(n_*maxValue_*(1.0 - sqr(coord)));
    vectorField::operator=(n_*0.374/0.4187*(log(coord)));
}


// Write
void functionVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("maxValue")
        << maxValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, functionVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
