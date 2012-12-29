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
    uStar_(0),
    kappa_(0),
    yZero_(0),
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
    uStar_(ptf.uStar_),
    kappa_(ptf.kappa_),
    yZero_(ptf.yZero_),
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
    uStar_(readScalar(dict.lookup("uStar"))),
    kappa_(readScalar(dict.lookup("kappa"))),
    yZero_(readScalar(dict.lookup("yZero"))),
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
    uStar_(fcvpvf.uStar_),
    kappa_(fcvpvf.kappa_),
    yZero_(fcvpvf.yZero_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void functionVelocityFvPatchVectorField::updateCoeffs()
{

//    Info << "Function called." << endl;

    // Get range and orientation
    //boundBox bb(patch().patch().localPoints(), false);
    //vector ctr = 0.5*(bb.max() + bb.min());

    const vectorField& c = patch().Cf();

    // Calculate argument for logarithmic profile
    scalarField coord = ((yZero_ + (c & y_))/yZero_);
    
    vectorField::operator=(n_*uStar_/kappa_*(log(coord)));
}


// Write
void functionVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("uStar")
        << uStar_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("yZero")
        << yZero_ << token::END_STATEMENT << nl;
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
