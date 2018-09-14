/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "packingLimitModel.H"
#include "SortableList.H"
#include "zeroGradientFvPatchFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
    defineTypeNameAndDebug(packingLimitModel, 0);

    defineRunTimeSelectionTable(packingLimitModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModel::packingLimitModel
(
    const dictionary& dict,
    const multiphaseKineticTheorySystem& kt
)
:
    constantDiameters_(dict.lookupOrDefault("constantDiameters", true)),
    dict_(dict),
    kt_(kt),
    mesh_(kt.fluid().mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::packingLimitModel::~packingLimitModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::packingLimitModel::alphaMax() const
{
    tmp<volScalarField> tmpMaxAlpha
    (
        new volScalarField
        (
            IOobject
            (
                "maxAlpha",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("0", dimless, 0.0),
            wordList
            (
                mesh_.boundaryMesh().size(),
                zeroGradientFvPatchScalarField::typeName
            )
        )
    );
    volScalarField& maxAlpha(tmpMaxAlpha.ref());

    const wordList& phases = kt_.phases();

    if (constantDiameters_)
    {
        // Sort diameters from largest to smallest
        scalarList ds(phases.size());
        forAll(phases, phasei)
        {
            ds[phasei] = kt_.fluid().phases()[phases[phasei]].d()()[0];
        }

        forAll(maxAlpha, celli)
        {
            maxAlpha[celli] = alphaMax(celli, ds);
        }
    }
    else
    {
        PtrList<volScalarField> dList(phases.size());
        forAll(phases, phasei)
        {
            dList.set
            (
                phasei,
                new volScalarField(kt_.fluid().phases()[phases[phasei]])
            );
        }

        forAll(maxAlpha, celli)
        {
            // Sort diameters from largest to smallest
            scalarList ds(phases.size());

            forAll(phases, phasei)
            {
                ds[phasei] = dList[phasei][celli];
            }

            maxAlpha[celli] = alphaMax(celli, ds);
        }
    }

    return tmpMaxAlpha;
}
// ************************************************************************* //
