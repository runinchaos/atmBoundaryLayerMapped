/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "atmBoundaryLayerMappedKFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerMappedKFvPatchScalarField::
atmBoundaryLayerMappedKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayerMapped(iF.time(), p.patch())
{}


atmBoundaryLayerMappedKFvPatchScalarField::
atmBoundaryLayerMappedKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayerMapped(iF.time(), p.patch(), dict, "k")
{
    phiName_ = dict.getOrDefault<word>("phi", "phi");

    refValue() = k(patch().Cf());
    refGrad() = 0;
    valueFraction() = 1;

    if (!initABL_)
    {
        scalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        scalarField::operator=(refValue());
        initABL_ = false;
    }
}


atmBoundaryLayerMappedKFvPatchScalarField::
atmBoundaryLayerMappedKFvPatchScalarField
(
    const atmBoundaryLayerMappedKFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(psf, p, iF, mapper),
    atmBoundaryLayerMapped(psf, p, mapper)
{}


atmBoundaryLayerMappedKFvPatchScalarField::
atmBoundaryLayerMappedKFvPatchScalarField
(
    const atmBoundaryLayerMappedKFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(psf, iF),
    atmBoundaryLayerMapped(psf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerMappedKFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Strict check: U field must exist and use atmBoundaryLayerMappedVelocity
    if (!db().foundObject<volVectorField>("U"))
    {
        FatalErrorInFunction
            << "U field not found in database. "
            << "atmBoundaryLayerMappedK must be used together with atmBoundaryLayerMappedVelocity."
            << abort(FatalError);
    }

    const fvPatchVectorField& Upatch = 
        db().lookupObject<volVectorField>("U").boundaryField()[patch().index()];

    if (Upatch.type() != "atmBoundaryLayerMappedVelocity")
    {
        FatalErrorInFunction
            << "U boundary field type '" << Upatch.type() << "' is not compatible. "
            << "atmBoundaryLayerMappedK must be used together with atmBoundaryLayerMappedVelocity, "
            << "but found '" << Upatch.type() << "' at patch '" << patch().name() << "'."
            << abort(FatalError);
    }

    // Get actual U values from patch
    const vectorField& Uvalues = Upatch.patchInternalField();

    // Calculate u* from actual U using base class function
    tmp<scalarField> tuStar = UstarFromU(Uvalues, patch().Cf());

    // Calculate k using the actual u* via base class function
    refValue() = kFromUstar(tuStar(), patch().Cf());

    inletOutletFvPatchScalarField::updateCoeffs();
}


void atmBoundaryLayerMappedKFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);
    atmBoundaryLayerMapped::autoMap(m);
}


void atmBoundaryLayerMappedKFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(psf, addr);

    const atmBoundaryLayerMappedKFvPatchScalarField& blpsf =
        refCast<const atmBoundaryLayerMappedKFvPatchScalarField>(psf);

    atmBoundaryLayerMapped::rmap(blpsf, addr);
}


void atmBoundaryLayerMappedKFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    atmBoundaryLayerMapped::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    atmBoundaryLayerMappedKFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
