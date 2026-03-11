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

#include "atmBoundaryLayerMappedVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerMappedVelocityFvPatchVectorField::
atmBoundaryLayerMappedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(p, iF),
    atmBoundaryLayerMapped(iF.time(), p.patch())
{}


atmBoundaryLayerMappedVelocityFvPatchVectorField::
atmBoundaryLayerMappedVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchVectorField(p, iF),
    atmBoundaryLayerMapped(iF.time(), p.patch(), dict, "U")
{
    phiName_ = dict.getOrDefault<word>("phi", "phi");

    refValue() = Umapped(patch().Cf());
    refGrad() = Zero;
    valueFraction() = 1;

    if (!initABL_)
    {
        vectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        vectorField::operator=(refValue());
        initABL_ = false;
    }
}


atmBoundaryLayerMappedVelocityFvPatchVectorField::
atmBoundaryLayerMappedVelocityFvPatchVectorField
(
    const atmBoundaryLayerMappedVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchVectorField(pvf, p, iF, mapper),
    atmBoundaryLayerMapped(pvf, p, mapper)
{}


atmBoundaryLayerMappedVelocityFvPatchVectorField::
atmBoundaryLayerMappedVelocityFvPatchVectorField
(
    const atmBoundaryLayerMappedVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    inletOutletFvPatchVectorField(pvf, iF),
    atmBoundaryLayerMapped(pvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerMappedVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    refValue() = Umapped(patch().Cf());

    inletOutletFvPatchVectorField::updateCoeffs();
}


void atmBoundaryLayerMappedVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchVectorField::autoMap(m);
    atmBoundaryLayerMapped::autoMap(m);
}


void atmBoundaryLayerMappedVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    inletOutletFvPatchVectorField::rmap(pvf, addr);

    const atmBoundaryLayerMappedVelocityFvPatchVectorField& blpvf =
        refCast<const atmBoundaryLayerMappedVelocityFvPatchVectorField>(pvf);

    atmBoundaryLayerMapped::rmap(blpvf, addr);
}


void atmBoundaryLayerMappedVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    atmBoundaryLayerMapped::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    atmBoundaryLayerMappedVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
