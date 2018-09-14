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

#include "multiphaseKineticTheorySystem.H"
#include "kineticTheoryModel.H"
#include "multiphaseSystem.H"
#include "packingLimitModel.H"
#include "radialModel.H"
#include "viscosityModel.H"
#include "frictionalStressModel.H"
#include "granularPressureModel.H"
#include "conductivityModel.H"
#include "phaseSystem.H"
#include "mathematicalConstants.H"
#include "SortableList.H"
#include "zeroGradientFvPatchFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseKineticTheorySystem::multiphaseKineticTheorySystem
(
    const phaseSystem& fluid
)
:
    regIOobject
    (
        IOobject
        (
            "kineticTheorySystem",
            fluid.mesh().time().constant(),
            fluid.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            true
        )
    ),
    fluid_(fluid),
    dict_(fluid.subDict("kineticTheory")),
    name_(dict_.lookup("name")),
    writeTotal_(dict_.lookupOrDefault("writeTotal", false)),
    alphap_
    (
        IOobject
        (
            IOobject::groupName("alpha", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimless, 0.0)
    ),
    Up_
    (
        IOobject
        (
            IOobject::groupName("U", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedVector("0", dimVelocity, Zero)
    ),
    Thetap_
    (
        IOobject
        (
            IOobject::groupName("Theta", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("0", sqr(dimVelocity), 0.0)
    ),
    phases_(),
    packingLimitModel_
    (
        kineticTheoryModels::packingLimitModel::New
        (
            dict_,
         *this
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            dict_,
            *this
        )
    ),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            dict_,
            *this
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
        (
            dict_,
            *this
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            dict_,
         *this
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            dict_,
            *this
        )
    ),
    eTable_
    (
        dict_.lookup("coeffRest")
    ),
    CfTable_
    (
        dict_.lookup("coeffFric")
    ),
    alphaMax_
    (
        IOobject
        (
            IOobject::groupName("alphaMax", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("one", dimless, 1.0),
        wordList
        (
            alphap_.boundaryField().size(),
            zeroGradientFvPatchScalarField::typeName
        )
    ),
    minAlphaMax_(1.0),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseKineticTheorySystem::~multiphaseKineticTheorySystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiphaseKineticTheorySystem::read()
{
    Info<<"multiphase read"<<endl;
    residualAlpha_.readIfPresent(dict_);

    radialModel_->read();
    viscosityModel_->read();
    frictionalStressModel_->read();
    granularPressureModel_->read();
    conductivityModel_->read();

    writeTotal_ = dict_.lookupOrDefault("writeTotal", false);
    Info<<"write Total: " << writeTotal_<<endl;
    if (writeTotal_)
    {
        alphap_.writeOpt() = AUTO_WRITE;
        Up_.writeOpt() = AUTO_WRITE;
        Thetap_.writeOpt() = AUTO_WRITE;
    }
    else
    {
        alphap_.writeOpt() = NO_WRITE;
        Up_.writeOpt() = NO_WRITE;
        Thetap_.writeOpt() = NO_WRITE;
    }

    return true;
}

bool Foam::multiphaseKineticTheorySystem::readIfModified()
{
    if (fluid_.modified())
    {
        read();

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::multiphaseKineticTheorySystem::polydisperse() const
{
    return (phases_.size() > 1);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseKineticTheorySystem::gs0
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return tmp<volScalarField>
    (
        *gs0Table_[phasePairKey(phase1.name(), phase2.name(), false)]
    );
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseKineticTheorySystem::gs0Prime
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return radialModel_->g0prime(phase1, phase2);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseKineticTheorySystem::nu
(
    const phaseModel& phase,
    const volScalarField& Theta
) const
{
    phasePairKey key(phase.name(), phase.name(), false);
    return viscosityModel_->nu
    (
        phase,
        Theta,
        *gs0Table_[key],
        phase.rho(),
        phase.d(),
        dimensionedScalar("e", dimless, eTable_[key])
    );
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseKineticTheorySystem::PsCoeff(const phaseModel& phase) const
{
    if (phases_.size() == 1)
    {
        phasePairKey key(phase.name(), phase.name(), false);
        const volScalarField& Theta =
            fluid_.mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("Theta", phase.name())
            );

        return tmp<volScalarField>
        (
            *PsCoeffs_[key]
           /max
            (
                Theta,
                dimensionedScalar("SMALL", sqr(dimVelocity), 1e-6)
            )
        );
    }

    tmp<volScalarField> tmpPsCoeff
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("PsCoeff", phase.name()),
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar("0", dimDensity, 0.0)
        )
    );
    volScalarField& psCoeff = tmpPsCoeff.ref();

    forAll(phases_, phasei)
    {
        const phaseModel& phase2 = fluid_.phases()[phases_[phasei]];
        phasePairKey key(phase.name(), phase2.name(), false);
        const volScalarField& Theta =
            fluid_.mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("Theta", phase.name())
            );

        psCoeff +=
            *PsCoeffs_[key]
           /max(Theta, dimensionedScalar("SMALL", sqr(dimVelocity), 1e-6));
    }
    return tmpPsCoeff;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseKineticTheorySystem::PsCoeffPrime(const phaseModel& phase) const
{
    if (phases_.size() == 1)
    {
        phasePairKey key(phase.name(), phase.name(), false);
        const volScalarField& Theta =
        fluid_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("Theta", phase.name())
        );

        return tmp<volScalarField>
        (
            *PsCoeffsPrime_[key]
           /max
            (
                Theta,
                dimensionedScalar("SMALL", sqr(dimVelocity), 1e-6)
            )
        );
    }

    tmp<volScalarField> tmpPsCoeffPrime
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("PsCoeffPrime", phase.name()),
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar("0", dimDensity, 0.0)
        )
    );
    volScalarField& psCoeffPrime = tmpPsCoeffPrime.ref();

    forAll(phases_, phasei)
    {
        const phaseModel& phase2 = fluid_.phases()[phases_[phasei]];
        phasePairKey key(phase.name(), phase2.name(), false);
        const volScalarField& Theta =
            fluid_.mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("Theta", phase.name())
            );

        psCoeffPrime +=
            *PsCoeffsPrime_[key]
           /max(Theta, dimensionedScalar("SMALL", sqr(dimVelocity), 1e-6));
    }
    return tmpPsCoeffPrime;
}

Foam::tmp<Foam::volScalarField>
Foam::multiphaseKineticTheorySystem::kappa
(
    const phaseModel& phase,
    const volScalarField& Theta
) const
{
    phasePairKey key(phase.name(), phase.name(), false);
    return conductivityModel_->kappa
    (
        phase,
        Theta,
        *gs0Table_[key],
        phase.rho(),
        phase.d(),
        dimensionedScalar("e", dimless, eTable_[key])
    );
}

Foam::tmp<Foam::volScalarField>
Foam::multiphaseKineticTheorySystem::
frictionalPressure(const phaseModel& phase) const
{
    return frictionalStressModel_->frictionalPressure
    (
        phase,
        alphaMax_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseKineticTheorySystem::
frictionalPressurePrime(const phaseModel& phase) const
{
    return frictionalStressModel_->frictionalPressurePrime
    (
        phase,
        alphaMax_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseKineticTheorySystem::nuFrictional(const phaseModel& phase) const
{
    tmp<volTensorField> tgradU(fvc::grad(Up_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    return frictionalStressModel_->nu
    (
        phase,
        alphaMax_,
        frictionalPressure(phase)/phase.rho(),
        D
    );
}



const Foam::wordList& Foam::multiphaseKineticTheorySystem::phases() const
{
    return phases_;
}


void Foam::multiphaseKineticTheorySystem::addPhase
(
    const phaseModel& phase
)
{
    const word& phaseName(phase.name());
    phases_.append(phaseName);

    minAlphaMax_ = min(minAlphaMax_, phase.alphaMax());

    // Print granular quantities only if more than 1 phase is present
    if (phases_.size() > 1 && writeTotal_)
    {
        alphap_.writeOpt() = AUTO_WRITE;
        Up_.writeOpt() = AUTO_WRITE;
        Thetap_.writeOpt() = AUTO_WRITE;
    }

    forAll(phases_, phasei)
    {
        phasePairKey key
        (
            phaseName,
            phases_[phasei],
            false
        );
        pairs_.append(key);
        word name = key.first() + "And" + key.second();

        if (phaseName == phases_[phasei])
        {
            name = phaseName;
        }

        gs0Table_.insert
        (
            key,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "gs0",
                        name
                    ),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedScalar("gs0", dimless, 0.0)
            )
        );

        PsCoeffs_.insert
        (
            key,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "PsCoeff",
                        name
                    ),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedScalar("PsCoeff", dimDensity*sqr(dimVelocity), 0.0)
            )
        );

        PsCoeffsPrime_.insert
        (
            key,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "PsCoeffPrime",
                        name
                    ),
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh()
                ),
                fluid_.mesh(),
                dimensionedScalar
                (
                    "PsCoeffPrime",
                    dimDensity*sqr(dimVelocity),
                    0.0
                )
            )
        );
    }
}


bool Foam::multiphaseKineticTheorySystem::found(const word& phaseName) const
{
    forAll(phases_, phasei)
    {
        if (phases_[phasei] == phaseName)
        {
            return true;
        }
    }
    return false;
}


void Foam::multiphaseKineticTheorySystem::correct()
{
    alphap_ = 0.0;
    Up_ = dimensionedVector("0", dimVelocity, Zero);
    Thetap_ = dimensionedScalar("0", sqr(dimVelocity), 0.0);

    forAll(phases_, phasei)
    {
        const phaseModel& phase = fluid_.phases()[phases_[phasei]];
        const volScalarField& alpha = fluid_.phases()[phases_[phasei]];
        const volScalarField& Theta =
            fluid_.mesh().lookupObject<volScalarField>
            (
                IOobject::groupName
                (
                    "Theta",
                    phase.name()
                )
            );
        alphap_ += alpha;
        Up_ += alpha*phase.U();
        Thetap_ +=  alpha*Theta;
    }
    Up_ /= max(alphap_, residualAlpha_);
    Thetap_ /= max(alphap_, residualAlpha_);

    forAll(pairs_, pairi)
    {
        const phasePairKey& key = pairs_[pairi];
        const phaseModel& phase1 = fluid_.phases()[key.first()];
        const phaseModel& phase2 = fluid_.phases()[key.second()];
        const volScalarField& Theta1 =
            fluid_.mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("Theta", phase1.name())
            );
        const volScalarField& Theta2 =
            fluid_.mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("Theta", phase2.name())
            );

        *gs0Table_[key] = radialModel_->g0(phase1, phase2);
        *PsCoeffs_[key] =
            granularPressureModel_->granularPressureCoeff
            (
                phase1,
                phase2,
                Theta1,
                Theta2,
                gs0(phase1, phase2),
                eTable_[key]
            );

        *PsCoeffsPrime_[key] =
            granularPressureModel_->granularPressureCoeffPrime
            (
                phase1,
                phase2,
                Theta1,
                Theta2,
                gs0(phase1, phase2),
                gs0Prime(phase1, phase2),
                eTable_[key]
            );
    }

    alphaMax_ = packingLimitModel_->alphaMax();
    alphaMax_.correctBoundaryConditions();
}


void Foam::multiphaseKineticTheorySystem::correctAlphap()
{
    alphap_ = 0.0;
    forAll(phases_, phasei)
    {
        alphap_ += fluid_.phases()[phases_[phasei]];
    }
    alphap_.correctBoundaryConditions();
}
// ************************************************************************* //
