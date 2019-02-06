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

#include "Syamlal.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "kineticTheorySystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Syamlal, 0);
    addToRunTimeSelectionTable(dragModel, Syamlal, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Syamlal::Syamlal
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    C1_
    (
        dimensionedScalar::lookupOrDefault
        (
            "C1",
            dict,
            dimensionSet(0, -2, 1, 0, 0, 0, 0),
            0.3
        )
    ),
    kineticTheorySystem_
    (
        pair_.phase1().mesh().lookupObject<kineticTheorySystem>
        (
            "kineticTheorySystem"
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Syamlal::~Syamlal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::Syamlal::CdRe() const
{
    FatalErrorInFunction
        << "Not implemented."
        << "Drag coefficient not defined for the Syamlal model."
        << exit(FatalError);

    return pair_.phase1();
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::Syamlal::K() const
{
    const phaseModel& phase1 = pair_.phase1();
    const phaseModel& phase2 = pair_.phase2();

    const volScalarField& Pfric = kineticTheorySystem_.frictionalPressure();
    scalar e = kineticTheorySystem_.es()[pair_];
    scalar Cf = kineticTheorySystem_.Cfs()[pair_];
    scalar pi = Foam::constant::mathematical::pi;

    tmp<volScalarField> g0 = kineticTheorySystem_.gs0(phase1, phase2);
    return
        (
            3.0*(1.0 + e)*(pi/2.0 + Cf*sqr(pi)/8.0)
           *phase1*phase1.rho()*phase2*phase2.rho()
           *sqr(phase1.d() + phase2.d())*g0*pair_.magUr()
        )
       /(
            2.0*pi
           *(
                phase1.rho()*pow3(phase1.d())
              + phase2.rho()*pow3(phase2.d())
            )
        )
      + C1_*Pfric;
}


Foam::tmp<Foam::surfaceScalarField> Foam::dragModels::Syamlal::Kf() const
{
    return fvc::interpolate(K());
}
// ************************************************************************* //
