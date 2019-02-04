/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "SchneiderbauerFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
    defineTypeNameAndDebug(Schneiderbauer, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        Schneiderbauer,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Schneiderbauer::Schneiderbauer
(
    const dictionary& dict
)
:
    frictionalStressModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    phi1_("phi1", dimless, coeffDict_),
    phi2_("phi2", dimless, coeffDict_),
    Io_("Io", dimless, coeffDict_),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        coeffDict_
    )
{
    phi1_ *= constant::mathematical::pi/180.0;
    phi2_ *= constant::mathematical::pi/180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Schneiderbauer::~Schneiderbauer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schneiderbauer::
frictionalPressure
(
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{

    return 4.0*rho*sqr(b_*ds*tr(dev(D))/(alphaMax - alphap));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schneiderbauer::
frictionalPressurePrime
(
    const volScalarField& alphap,
    const volScalarField& alphaMax
) const
{
    volScalarField alphaMinFriction(alphaMinFrictionByAlphap_*alphaMax);

    return 8.0*rho*sqr(b_*ds*tr(dev(D)))/pow3(alphaMax - alphap);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schneiderbauer::nu
(
    const phaseModel& phase,
    const volScalarField& alphap,
    const volScalarField& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D
) const
{
    volScalarField alphaMinFriction(alphaMinFrictionByAlphap_*alphaMax);

    tmp<volScalarField> tnu
    (
        new volScalarField
        (
            IOobject
            (
                "Schneiderbauer:nu",
                phase.mesh().time().timeName(),
                phase.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phase.mesh(),
            dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
        )
    );

    volScalarField& nuf = tnu.ref();
    volScalarField Is(2.0*tr(D)*ds/sqrt(pf/rho));
    volScalarField mui(tan(phi1_) + (tan(phi2_) - tan(phi1_))/(Io_/Is + 1.0));


    forAll(D, celli)
    {
        if (alphap[celli] > alphaMinFriction[celli])
        {
            nuf[celli] = mui[celli]*pf[celli]/(2.0*tr(D));
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    const volVectorField& U = phase.U();

    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            nufBf[patchi] =
                (
                    pf.boundaryField()[patchi]*sin(phi_.value())
                   /(
                        mag(U.boundaryField()[patchi].snGrad())
                      + SMALL
                    )
                );
        }
    }

    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
}


bool Foam::kineticTheoryModels::frictionalStressModels::Schneiderbauer::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    phi_.read(coeffDict_);
    phi_ *= constant::mathematical::pi/180.0;
    alphaMinFriction_.read(coeffDict_);

    return true;
}


// ************************************************************************* //
