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

#include "HuilinPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
    defineTypeNameAndDebug(Huilin, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        Huilin,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Huilin::Huilin
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    granularPressureModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Huilin::~Huilin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Huilin::granularPressureCoeff
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& Theta1,
    const volScalarField& Theta2,
    const volScalarField& g0,
    const dimensionedScalar& e
) const
{
    const scalar pi = Foam::constant::mathematical::pi;
    const volScalarField& alpha1 = phase1;
    const volScalarField& alpha2 = phase2;
    const volScalarField& rho1 = phase1.rho();
    const volScalarField& rho2 = phase2.rho();
    volScalarField d1(phase1.d());
    volScalarField d2(phase2.d());
    volScalarField d12(0.5*(d1 + d2));

    volScalarField m1(pi/6.0*pow3(d1)*phase1.rho());
    volScalarField m2(pi/6.0*pow3(d2)*phase2.rho());
    volScalarField n1(6.0*phase1/(pi*pow3(d1)));
    volScalarField n2(6.0*phase2/(pi*pow3(d2)));
    volScalarField m0(m1 + m2);
    volScalarField omega
    (
        (m1*Theta1 - m2*Theta2)
       /sqrt
        (
            sqr(m1)*sqr(Theta1) + sqr(m2)*sqr(Theta2)
          + Theta1*Theta2*(sqr(m1) + sqr(m2))
          + dimensionedScalar("0", sqr(dimMass*sqr(dimVelocity)), 1e-10)
        )
    );

    return
    (
        pi*(1 + e)*pow3(d12)*g0*alpha1*rho1*alpha2*rho2*m0*Theta1*Theta2
       /(
            3.0*(sqr(m1)*Theta1 + sqr(m2)*Theta2)
          + dimensionedScalar("0", sqr(dimMass*dimVelocity), 1e-10)
        )
       *pow
        (
            sqr(m0)*Theta1*Theta2
           /(
                (sqr(m1)*Theta1 + sqr(m2)*Theta2)*(Theta1 + Theta2)
              + dimensionedScalar("0", sqr(dimMass*sqr(dimVelocity)), 1e-10)
            ),
            3.0/2.0
        )
       *(1.0 - 3.0*omega + 6.0*sqr(omega))
    );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Huilin::
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
    const scalar pi = Foam::constant::mathematical::pi;
    volScalarField d1(phase1.d());
    volScalarField d2(phase2.d());
    volScalarField d12(0.5*(d1 + d2));

    volScalarField m1(pi/6.0*pow3(d1)*phase1.rho());
    volScalarField m2(pi/6.0*pow3(d2)*phase2.rho());
    volScalarField n1(6.0*phase1/(pi*pow3(d1)));
    volScalarField n2(6.0*phase2/(pi*pow3(d2)));
    volScalarField m0(m1 + m2);
    volScalarField omega
    (
        (m1*Theta1 - m2*Theta2)
       /sqrt
        (
            sqr(m1)*sqr(Theta1) + sqr(m2)*sqr(Theta2)
          + Theta1*Theta2*(sqr(m1) + sqr(m2))
          + dimensionedScalar("0", sqr(dimMass*sqr(dimVelocity)), 1e-10)
        )
    );

    volScalarField granularPressurePrime
    (
        pi*(1 + e)*pow3(d12)*n2*m2*m0*Theta1*Theta2
       /(
            3.0*(sqr(m1)*Theta1 + sqr(m2)*Theta2)
          + dimensionedScalar("0", sqr(dimMass*dimVelocity), 1e-10)
        )
       *pow
        (
            sqr(m0)*Theta1*Theta2
           /(
                (sqr(m1)*Theta1 + sqr(m2)*Theta2)*(Theta1 + Theta2)
              + dimensionedScalar("0", sqr(dimMass*sqr(dimVelocity)), 1e-10)
            ),
            3.0/2.0
        )
       *(1.0 - 3.0*omega + 6.0*sqr(omega))
    );

    return granularPressurePrime*(g0prime*m1*n1 + g0*phase1.rho());
}


// ************************************************************************* //
