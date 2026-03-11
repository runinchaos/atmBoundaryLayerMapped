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

#include "atmBoundaryLayerMappedEpsilonFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerMappedEpsilonFvPatchScalarField::
atmBoundaryLayerMappedEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayerMapped(iF.time(), p.patch())
{}


atmBoundaryLayerMappedEpsilonFvPatchScalarField::
atmBoundaryLayerMappedEpsilonFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayerMapped(iF.time(), p.patch(), dict)
{
    phiName_ = dict.getOrDefault<word>("phi", "phi");

    refValue() = epsilon(patch().Cf());
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


atmBoundaryLayerMappedEpsilonFvPatchScalarField::
atmBoundaryLayerMappedEpsilonFvPatchScalarField
(
    const atmBoundaryLayerMappedEpsilonFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(psf, p, iF, mapper),
    atmBoundaryLayerMapped(psf, p, mapper)
{}


atmBoundaryLayerMappedEpsilonFvPatchScalarField::
atmBoundaryLayerMappedEpsilonFvPatchScalarField
(
    const atmBoundaryLayerMappedEpsilonFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(psf, iF),
    atmBoundaryLayerMapped(psf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerMappedEpsilonFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Check if U field exists and is using mapped velocity BC
    const volVectorField* Uptr = nullptr;
    bool useMappedU = false;

    if (db().foundObject<volVectorField>("U"))
    {
        Uptr = &db().lookupObject<volVectorField>("U");
        const fvPatchVectorField& Upatch = Uptr->boundaryField()[patch().index()];

        // Check if U is using atmBoundaryLayerMappedVelocity
        if (Upatch.type() == "atmBoundaryLayerMappedVelocity")
        {
            useMappedU = true;
        }
    }

    if (useMappedU && Uptr)
    {
        // Get actual U values from patch
        const fvPatchVectorField& Upatch = Uptr->boundaryField()[patch().index()];
        const vectorField& Uvalues = Upatch.patchInternalField();

        // Calculate u* from actual U and then epsilon
        tmp<scalarField> tuStar = UstarFromU(Uvalues, patch().Cf());
        refValue() = epsilonFromUstar(tuStar(), patch().Cf());
    }
    else
    {
        // Use standard ABL formula based on Uref
        refValue() = epsilon(patch().Cf());
    }

    inletOutletFvPatchScalarField::updateCoeffs();
}


void atmBoundaryLayerMappedEpsilonFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);
    atmBoundaryLayerMapped::autoMap(m);
}


void atmBoundaryLayerMappedEpsilonFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(psf, addr);

    const atmBoundaryLayerMappedEpsilonFvPatchScalarField& blpsf =
        refCast<const atmBoundaryLayerMappedEpsilonFvPatchScalarField>(psf);

    atmBoundaryLayerMapped::rmap(blpsf, addr);
}


void atmBoundaryLayerMappedEpsilonFvPatchScalarField::write(Ostream& os) const
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
    atmBoundaryLayerMappedEpsilonFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
