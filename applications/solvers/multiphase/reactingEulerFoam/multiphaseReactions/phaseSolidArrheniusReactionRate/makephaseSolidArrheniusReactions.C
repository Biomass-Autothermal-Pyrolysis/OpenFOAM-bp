/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "makeReaction.H"
#include "reactionTypes.H"
#include "phaseSolidArrheniusReactionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeGeneralReaction
    (
        gasHThermoPhysics,
        IrreversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        gasHThermoPhysics,
        ReversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        gasEThermoPhysics,
        IrreversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        gasEThermoPhysics,
        ReversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        constGasHThermoPhysics,
        IrreversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        constGasHThermoPhysics,
        ReversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        constGasEThermoPhysics,
        IrreversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        constGasEThermoPhysics,
        ReversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        constHThermoPhysics,
        IrreversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        constHThermoPhysics,
        ReversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        constEThermoPhysics,
        IrreversibleReaction,
        phaseSolidArrheniusReactionRate
    )

    makeGeneralReaction
    (
        constEThermoPhysics,
        ReversibleReaction,
        phaseSolidArrheniusReactionRate
    )
}

// ************************************************************************* //
