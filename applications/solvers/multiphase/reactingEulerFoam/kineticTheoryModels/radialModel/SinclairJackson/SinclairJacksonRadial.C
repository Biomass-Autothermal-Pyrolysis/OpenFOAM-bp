/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "SinclairJacksonRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(SinclairJackson, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        SinclairJackson,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::SinclairJackson::SinclairJackson
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    radialModel(dict, kt),
    alphaMinFriction_(dict.lookupType<scalar>("alphaMinFriction"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::SinclairJackson::~SinclairJackson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::SinclairJackson::g0
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return
        1.0
       /(
           1.0
         - cbrt(min(phase1, alphaMinFriction_)/phase1.alphaMax())
        );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::SinclairJackson::g0prime
(
    const phaseModel& phase1,
 const phaseModel& phase2
) const
{
    volScalarField aByaMax
    (
        cbrt
        (
            min
            (
                max(phase1, scalar(1e-3)),
                alphaMinFriction_
            )/phase1.alphaMax()
        )
    );

    return (1.0/(3*phase1.alphaMax()))/sqr(aByaMax - sqr(aByaMax));
}


// ************************************************************************* //
