/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "atmBoundaryLayerMapped.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerMapped::atmBoundaryLayerMapped(const Time& time, const polyPatch& pp)
:
    initABL_(false),
    kappa_(0.41),
    Cmu_(0.09),
    C1_(0.0),
    C2_(1.0),
    ppMin_((boundBox(pp.points())).min()),
    time_(time),
    patch_(pp),
    zDir_(nullptr),
    Zref_(nullptr),
    z0_(nullptr),
    d_(nullptr),
    UMapper_(nullptr),
    scalarMapper_(nullptr)
{}


atmBoundaryLayerMapped::atmBoundaryLayerMapped
(
    const Time& time,
    const polyPatch& pp,
    const dictionary& dict,
    const word& fieldName
)
:
    initABL_(dict.getOrDefault<bool>("initABL", true)),
    kappa_
    (
        dict.getCheckOrDefault<scalar>("kappa", 0.41, scalarMinMax::ge(SMALL))
    ),
    Cmu_(dict.getCheckOrDefault<scalar>("Cmu", 0.09, scalarMinMax::ge(SMALL))),
    C1_(dict.getOrDefault("C1", 0.0)),
    C2_(dict.getOrDefault("C2", 1.0)),
    ppMin_((boundBox(pp.points())).min()),
    time_(time),
    patch_(pp),
    zDir_(Function1<vector>::New("zDir", dict, &time)),
    Zref_(Function1<scalar>::New("Zref", dict, &time)),
    z0_(PatchFunction1<scalar>::New(pp, "z0", dict)),
    d_(PatchFunction1<scalar>::New(pp, "d", dict)),
    UMapper_(nullptr),
    scalarMapper_(nullptr)
{
    // Always use mapping mode - this BC is designed for mapped data only
    // Force mapMethod to nearest to avoid triangulation issues
    const_cast<dictionary&>(dict).set("mapMethod", "nearest");
    if (fieldName == "U")
    {
        UMapper_.reset
        (
            new PatchFunction1Types::MappedFile<vector>
            (
                pp,
                "MappedFile",  // redirectType
                fieldName,     // entryName
                dict,          // dictionary
                true           // faceValues
            )
        );
    }
    else
    {
        scalarMapper_.reset
        (
            new PatchFunction1Types::MappedFile<scalar>
            (
                pp,
                "MappedFile",  // redirectType
                fieldName,     // entryName
                dict,          // dictionary
                true           // faceValues
            )
        );
    }
}


atmBoundaryLayerMapped::atmBoundaryLayerMapped
(
    const atmBoundaryLayerMapped& abl,
    const fvPatch& patch,
    const fvPatchFieldMapper& mapper
)
:
    initABL_(abl.initABL_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    C1_(abl.C1_),
    C2_(abl.C2_),
    ppMin_(abl.ppMin_),
    time_(abl.time_),
    patch_(patch.patch()),
    zDir_(abl.zDir_.clone()),
    Zref_(abl.Zref_.clone()),
    z0_(abl.z0_.clone(patch_)),
    d_(abl.d_.clone(patch_)),
    UMapper_(abl.UMapper_.clone(patch_)),
    scalarMapper_(abl.scalarMapper_.clone(patch_))
{}


atmBoundaryLayerMapped::atmBoundaryLayerMapped(const atmBoundaryLayerMapped& abl)
:
    initABL_(abl.initABL_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    C1_(abl.C1_),
    C2_(abl.C2_),
    ppMin_(abl.ppMin_),
    time_(abl.time_),
    patch_(abl.patch_),
    zDir_(abl.zDir_.clone()),
    Zref_(abl.Zref_.clone()),
    z0_(abl.z0_.clone(patch_)),
    d_(abl.d_.clone(patch_)),
    UMapper_(abl.UMapper_.clone(patch_)),
    scalarMapper_(abl.scalarMapper_.clone(patch_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector atmBoundaryLayerMapped::zDir() const
{
    const scalar t = time_.timeOutputValue();
    const vector dir(zDir_->value(t));
    const scalar magDir = mag(dir);

    if (magDir < SMALL)
    {
        FatalErrorInFunction
            << "magnitude of " << zDir_->name() << " = " << magDir
            << abort(FatalError);
    }

    return dir/magDir;
}


tmp<scalarField> atmBoundaryLayerMapped::z0() const
{
    const scalar t = time_.timeOutputValue();
    return max(z0_->value(t), ROOTVSMALL);
}


tmp<scalarField> atmBoundaryLayerMapped::d() const
{
    const scalar t = time_.timeOutputValue();
    return d_->value(t);
}


tmp<scalarField> atmBoundaryLayerMapped::UstarFromU
(
    const vectorField& Uvalues,
    const vectorField& pCf
) const
{
    // Safety checks
    if (!zDir_)
    {
        FatalErrorInFunction << "zDir not initialized" << abort(FatalError);
    }
    if (!z0_)
    {
        FatalErrorInFunction << "z0 not initialized" << abort(FatalError);
    }
    if (!d_)
    {
        FatalErrorInFunction << "d not initialized" << abort(FatalError);
    }
    if (Uvalues.size() != pCf.size())
    {
        FatalErrorInFunction << "Uvalues size " << Uvalues.size() << " != pCf size " << pCf.size() << abort(FatalError);
    }

    const scalar t = time_.timeOutputValue();
    const scalarField dvals(d_->value(t));
    const scalarField z0vals(max(z0_->value(t), ROOTVSMALL));
    const scalar groundMin = zDir() & ppMin_;

    // Calculate flow direction from Uvalues using parallel-safe reduce
    vector avgU = Zero;
    forAll(Uvalues, i)
    {
        avgU += Uvalues[i];
    }

    // Parallel: sum across all processors
    reduce(avgU, sumOp<vector>());

    // Get total number of faces across all processors
    label nTotalFaces = Uvalues.size();
    reduce(nTotalFaces, sumOp<label>());

    if (nTotalFaces > 0)
    {
        avgU /= nTotalFaces;
    }

    vector flowDirection;
    const scalar magAvgU = mag(avgU);
    if (magAvgU > SMALL)
    {
        flowDirection = avgU / magAvgU;
    }
    else
    {
        flowDirection = vector(1, 0, 0); // Default to x-direction
    }

    // Calculate u* from actual U using log law:
    // u* = kappa * U / ln((z - d + z0)/z0)
    scalarField ustar(Uvalues.size());

    forAll(Uvalues, i)
    {
        const scalar z = ((zDir() & pCf[i]) - groundMin);
        const scalar z0i = z0vals[i];
        const scalar di = dvals[i];

        // Get streamwise velocity component
        const scalar Ustream = Uvalues[i] & flowDirection;

        // Ensure log argument is valid (> 0)
        const scalar logArg = max((z - di + z0i)/z0i, ROOTVSMALL);

        if (Ustream > SMALL && logArg > 1.0)
        {
            ustar[i] = kappa_*Ustream/log(logArg);
        }
        else
        {
            ustar[i] = ROOTVSMALL;
        }
    }

    return ustar;
}


void atmBoundaryLayerMapped::autoMap(const fvPatchFieldMapper& mapper)
{
    if (z0_)
    {
        z0_->autoMap(mapper);
    }
    if (d_)
    {
        d_->autoMap(mapper);
    }
    if (UMapper_)
    {
        UMapper_->autoMap(mapper);
    }
    if (scalarMapper_)
    {
        scalarMapper_->autoMap(mapper);
    }
}


void atmBoundaryLayerMapped::rmap
(
    const atmBoundaryLayerMapped& abl,
    const labelList& addr
)
{
    if (z0_)
    {
        z0_->rmap(abl.z0_(), addr);
    }
    if (d_)
    {
        d_->rmap(abl.d_(), addr);
    }
    if (UMapper_)
    {
        UMapper_->rmap(abl.UMapper_(), addr);
    }
    if (scalarMapper_)
    {
        scalarMapper_->rmap(abl.scalarMapper_(), addr);
    }
}


// Note: U(), k(), epsilon(), omega() functions removed
// This BC only supports mapped mode, use Umapped() instead

tmp<scalarField> atmBoundaryLayerMapped::kFromUstar(const scalarField& uStar, const vectorField& pCf) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));
    const scalar groundMin = zDir() & ppMin_;

    // (YGCJ:Eq. 21) with bounds checking
    scalarField logArg = ((zDir() & pCf) - groundMin - d + z0)/z0;
    logArg = max(logArg, ROOTVSMALL);

    return sqr(uStar)/sqrt(Cmu_)
       *sqrt(C1_*log(logArg) + C2_);
}


tmp<scalarField> atmBoundaryLayerMapped::epsilonFromUstar(const scalarField& uStar, const vectorField& pCf) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));
    const scalar groundMin = zDir() & ppMin_;

    // (YGCJ:Eq. 22) with bounds checking
    scalarField denom = max((zDir() & pCf) - groundMin - d + z0, ROOTVSMALL);
    scalarField logArg = denom/z0;
    logArg = max(logArg, ROOTVSMALL);

    return pow3(uStar)/(kappa_*denom)
       *sqrt(C1_*log(logArg) + C2_);
}


tmp<scalarField> atmBoundaryLayerMapped::omegaFromUstar(const scalarField& uStar, const vectorField& pCf) const
{
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));
    const scalar groundMin = zDir() & ppMin_;

    // (YGJ:Eq. 13) with bounds checking
    scalarField denom = max((zDir() & pCf) - groundMin - d + z0, ROOTVSMALL);

    return uStar/(kappa_*sqrt(Cmu_)*denom);
}


tmp<vectorField> atmBoundaryLayerMapped::Umapped(const vectorField& pCf) const
{
    // Always use mapping - this BC is designed for mapped data only
    const scalar t = time_.timeOutputValue();
    tmp<vectorField> tmappedU(UMapper_->value(t));
    vectorField& mappedU = tmappedU.ref();

    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));
    const scalar groundMin = zDir() & ppMin_;
    const scalar ZrefVal = Zref_->value(t);

    const scalarField zHeight = (zDir() & pCf) - groundMin;

    forAll(mappedU, i)
    {
        const scalar z = zHeight[i];

        // Only apply atm correction below Zref; use mapped value above Zref
        if (z < ZrefVal)
        {
            const scalar z0i = z0[i];
            const scalar di = d[i];

            // Bounds checking for log arguments
            scalar logArg1 = max((z - di + z0i)/z0i, ROOTVSMALL);
            scalar logArg2 = max((ZrefVal + z0i)/z0i, ROOTVSMALL);

            scalar scale = log(logArg1) / log(logArg2);
            scale = max(scale, 0);

            mappedU[i] *= scale;
        }
    }

    return tmappedU;
}


tmp<scalarField> atmBoundaryLayerMapped::scalarMapped() const
{
    // Always use mapping - this BC is designed for mapped data only
    const scalar t = time_.timeOutputValue();
    return scalarMapper_->value(t);
}


void atmBoundaryLayerMapped::write(Ostream& os) const
{
    os.writeEntry("initABL", initABL_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("C1", C1_);
    os.writeEntry("C2", C2_);
    if (zDir_)
    {
        zDir_->writeData(os);
    }
    if (Zref_)
    {
        Zref_->writeData(os);
    }
    if (z0_)
    {
        z0_->writeData(os) ;
    }
    if (d_)
    {
        d_->writeData(os);
    }
    if (UMapper_)
    {
        UMapper_->writeData(os);
    }
    if (scalarMapper_)
    {
        scalarMapper_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
