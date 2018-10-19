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

#include "ChaoPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
    defineTypeNameAndDebug(Chao, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        Chao,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Chao::Chao
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    granularPressureModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Chao::~Chao()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Chao::granularPressureCoeff
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& Theta1,
    const volScalarField& Theta2,
    const volScalarField& g0,
    const dimensionedScalar& e
) const
{
    volScalarField m0
    (
        constant::mathematical::pi/6.0
       *(phase1.rho()*pow3(phase1.d()) + phase2.rho()*pow3(phase2.d()))
    );
    return
        phase1*phase2
       *constant::mathematical::pi/(3.0*m0)*phase1.rho()*phase2.rho()
       *(1.0 + e)*pow3((phase1.d() + phase2.d())/2.0)*g0
       *(Theta1 + Theta2 + 0.2*magSqr(phase1.U() - phase2.U()));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Chao::
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
    volScalarField m0
    (
        constant::mathematical::pi/6.0
       *(phase1.rho()*pow3(phase1.d()) + phase2.rho()*pow3(phase2.d()))
    );
    tmp<volScalarField> pCoeff
    (
        2.0*phase2*constant::mathematical::pi/(3.0*m0)*phase1.rho()
       *phase2.rho()*(1.0 + e)*pow3((phase1.d() + phase2.d())/2.0)
       *(Theta1 + Theta2 + 0.2*magSqr(phase1.U() - phase2.U()))
    );

    if (&phase1 == &phase2)
    {
        return pCoeff*(2.0*g0 + phase1*g0prime);

    }
    else
    {
        return pCoeff*(g0 + phase1*g0prime);
    }
}


// ************************************************************************* //
