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

#include "PrincetonPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
    defineTypeNameAndDebug(Princeton, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        Princeton,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Princeton::Princeton
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    granularPressureModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Princeton::~Princeton()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Princeton::granularPressureCoeff
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& Theta1,
    const volScalarField& Theta2,
    const volScalarField& g0,
    const dimensionedScalar& e
) const
{
    dimensionedScalar eta = 0.5*(1.0 + e);

    if (&phase1 == &phase2)
    {
        return phase1*phase1.rho()*Theta1
       *(
            1.0
          + 4.0*eta*phase1*g0
        );
    }
    return phase1*phase1.rho()*Theta1*4.0*eta*phase2*g0;
}



Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Princeton::
granularPressureCoeffPrime
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& Theta1,
    const volScalarField& Theta2,
    const volScalarField& g0,
    const volScalarField& g0prime,
    const dimensionedScalar& e
) const
{
    dimensionedScalar eta = 0.5*(1.0 + e);

    if (&phase1 == &phase2)
    {
        return phase1.rho()*Theta1
       *(
            1.0
          + 4.0*eta
           *(
               phase1*g0
             + phase1*phase2*g0prime
            )
        );
    }
    return
        phase1.rho()*Theta1*4.0*eta
       *(
           phase2*g0
         + phase1*phase2*g0prime
        );
}


// ************************************************************************* //
