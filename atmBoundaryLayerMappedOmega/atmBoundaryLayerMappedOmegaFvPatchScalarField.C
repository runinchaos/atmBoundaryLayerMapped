/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "atmBoundaryLayerMappedOmegaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerMappedOmegaFvPatchScalarField::
atmBoundaryLayerMappedOmegaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayerMapped(iF.time(), p.patch())
{}


atmBoundaryLayerMappedOmegaFvPatchScalarField::
atmBoundaryLayerMappedOmegaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayerMapped(iF.time(), p.patch(), dict, "omega")
{
    phiName_ = dict.getOrDefault<word>("phi", "phi");

    refValue() = 1;
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


atmBoundaryLayerMappedOmegaFvPatchScalarField::
atmBoundaryLayerMappedOmegaFvPatchScalarField
(
    const atmBoundaryLayerMappedOmegaFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(psf, p, iF, mapper),
    atmBoundaryLayerMapped(psf, p, mapper)
{}


atmBoundaryLayerMappedOmegaFvPatchScalarField::
atmBoundaryLayerMappedOmegaFvPatchScalarField
(
    const atmBoundaryLayerMappedOmegaFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(psf, iF),
    atmBoundaryLayerMapped(psf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerMappedOmegaFvPatchScalarField::updateCoeffs()
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
            << "atmBoundaryLayerMappedOmega must be used together with atmBoundaryLayerMappedVelocity."
            << abort(FatalError);
    }

    const fvPatchVectorField& Upatch =
        db().lookupObject<volVectorField>("U").boundaryField()[patch().index()];

    if (Upatch.type() != "atmBoundaryLayerMappedVelocity")
    {
        FatalErrorInFunction
            << "U boundary field type '" << Upatch.type() << "' is not compatible. "
            << "atmBoundaryLayerMappedOmega must be used together with atmBoundaryLayerMappedVelocity, "
            << "but found '" << Upatch.type() << "' at patch '" << patch().name() << "'."
            << abort(FatalError);
    }

    // Get actual U values from patch
    const vectorField Uvalues(Upatch.patchInternalField());

    // Calculate u* from actual U and then omega
    tmp<scalarField> tuStar = UstarFromU(Uvalues, patch().Cf());
    refValue() = omegaFromUstar(tuStar(), patch().Cf());

    inletOutletFvPatchScalarField::updateCoeffs();
}


void atmBoundaryLayerMappedOmegaFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);
    atmBoundaryLayerMapped::autoMap(m);
}


void atmBoundaryLayerMappedOmegaFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(psf, addr);

    const atmBoundaryLayerMappedOmegaFvPatchScalarField& blpsf =
        refCast<const atmBoundaryLayerMappedOmegaFvPatchScalarField>(psf);

    atmBoundaryLayerMapped::rmap(blpsf, addr);
}


void atmBoundaryLayerMappedOmegaFvPatchScalarField::write(Ostream& os) const
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
    atmBoundaryLayerMappedOmegaFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
