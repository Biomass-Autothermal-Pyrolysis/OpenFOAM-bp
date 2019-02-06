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

#include "kineticTheorySystem.H"
#include "kineticTheoryModel.H"
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

Foam::kineticTheorySystem::kineticTheorySystem
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
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        )
    ),
    fluid_(fluid),
    dict_(fluid.subDict("kineticTheory")),
    name_
    (
        dict_.lookupOrDefault<word>
        (
            "name",
            "kineticTheoryTotal"
        )
    ),
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
    kineticTheoryModels_(0),
    Thetas_(0),
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
        kineticTheoryModels::frictionalStressModel::New(dict_)
    ),
    eTable_(dict_.lookup("coeffRest")),
    CfTable_(dict_.lookup("coeffFric")),
    alphaMax_
    (
        IOobject
        (
            IOobject::groupName("alphaMax", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("one", dimless, 0.0),
        wordList
        (
            alphap_.boundaryField().size(),
            zeroGradientFvPatchScalarField::typeName
        )
    ),
    minAlphaMax_(1.0),
    residualAlpha_("residualAlpha", dimless, dict_),
    Pfr_
    (
        IOobject
        (
            IOobject::groupName("Pfr", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimPressure, 0.0)
    ),

    PfrPrime_
    (
        IOobject
        (
            IOobject::groupName("PfrPrime", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimPressure, 0.0)
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheorySystem::~kineticTheorySystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::kineticTheorySystem::read()
{
    residualAlpha_.readIfPresent(dict_);

    radialModel_->read();
    viscosityModel_->read();
    frictionalStressModel_->read();
    granularPressureModel_->read();
    conductivityModel_->read();

    writeTotal_ = dict_.lookupOrDefault("writeTotal", false);

    if (phaseIndexes_.size() > 1 && writeTotal_)
    {
        alphap_.writeOpt() = AUTO_WRITE;
        Up_.writeOpt() = AUTO_WRITE;
        Thetap_.writeOpt() = AUTO_WRITE;
        alphaMax_.writeOpt() = AUTO_WRITE;
    }
    else
    {
        alphap_.writeOpt() = NO_WRITE;
        Up_.writeOpt() = NO_WRITE;
        Thetap_.writeOpt() = NO_WRITE;
        alphaMax_.writeOpt() = NO_WRITE;
    }

    return true;
}

bool Foam::kineticTheorySystem::readIfModified()
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


bool Foam::kineticTheorySystem::polydisperse() const
{
    return (phaseIndexes_.size() > 1);
}


Foam::tmp<Foam::volScalarField> Foam::kineticTheorySystem::gs0
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


Foam::tmp<Foam::volScalarField> Foam::kineticTheorySystem::gs0Prime
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    return radialModel_->g0prime(phase1, phase2);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::nu
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
Foam::kineticTheorySystem::PsCoeff(const phaseModel& phase) const
{
    if (phaseIndexes_.size() == 1)
    {
        phasePairKey key(phase.name(), phase.name(), false);
        const volScalarField& Theta = Thetas_[phase.name()];

        return tmp<volScalarField>
        (
            *PsCoeffs_[key]
           /max
            (
                Theta,
                dimensionedScalar("SMALL", sqr(dimVelocity), 1e-8)
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

    forAll(phaseIndexes_, phasei)
    {
        const phaseModel& phase2 = fluid_.phases()[phaseIndexes_[phasei]];
        phasePairKey key(phase.name(), phase2.name(), false);
        const volScalarField& Theta = Thetas_[phase.name()];

        psCoeff +=
            *PsCoeffs_[key]
           /max(Theta, dimensionedScalar("SMALL", sqr(dimVelocity), 1e-8));
    }
    return tmpPsCoeff;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::PsCoeffPrime(const phaseModel& phase) const
{
    if (phaseIndexes_.size() == 1)
    {
        phasePairKey key(phase.name(), phase.name(), false);
        const volScalarField& Theta = Thetas_[phase.name()];

        return tmp<volScalarField>
        (
            *PsCoeffsPrime_[key]
           /max
            (
                Theta,
                dimensionedScalar("SMALL", sqr(dimVelocity), 1e-8)
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

    forAll(phaseIndexes_, phasei)
    {
        const phaseModel& phase2 = fluid_.phases()[phaseIndexes_[phasei]];
        phasePairKey key(phase.name(), phase2.name(), false);
        const volScalarField& Theta = Thetas_[phase.name()];

        psCoeffPrime +=
            *PsCoeffsPrime_[key]
           /max(Theta, dimensionedScalar("SMALL", sqr(dimVelocity), 1e-8));
    }
    return tmpPsCoeffPrime;
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheorySystem::kappa
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

const Foam::volScalarField&
Foam::kineticTheorySystem::frictionalPressure() const
{
    return Pfr_;
}


const Foam::volScalarField&
Foam::kineticTheorySystem::frictionalPressurePrime() const
{
    return PfrPrime_;
}


const Foam::volScalarField& Foam::kineticTheorySystem::nuFrictional() const
{
    return nuFric_;
}


const Foam::labelList& Foam::kineticTheorySystem::phaseIndexes() const
{
    return phaseIndexes_;
}


void Foam::kineticTheorySystem::addPhase
(
    RASModels::kineticTheoryModel& kt
)
{
    const phaseModel& phase = kt.phase();
    word phaseName(phase.name());
    kineticTheoryModels_.append(phaseName, new tmp<RASModels::kineticTheoryModel>(kt));
    Thetas_.append(phaseName, new tmp<volScalarField>(kt.Theta()));
    phaseIndexes_.append(phase.index());

    minAlphaMax_ = min(minAlphaMax_, phase.alphaMax());

    // Print granular quantities only if more than 1 phase is present
    if (phaseIndexes_.size() > 1 && writeTotal_)
    {
        alphap_.writeOpt() = AUTO_WRITE;
        Up_.writeOpt() = AUTO_WRITE;
        Thetap_.writeOpt() = AUTO_WRITE;
        alphaMax_.writeOpt() = AUTO_WRITE;
    }

    forAll(phaseIndexes_, phasei)
    {
        word otherPhaseName = kineticTheoryModels_[phasei]().phase().name();
        phasePairKey key
        (
            phaseName,
            otherPhaseName,
            false
        );
        pairs_.append(key);
        word name(key.second());
        name[0] = toupper(name[0]);
        name = key.first() + "And" + name;

        if (phaseName == otherPhaseName)
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


bool Foam::kineticTheorySystem::found(const word& phaseName) const
{
    forAll(phaseIndexes_, phasei)
    {
        if (fluid_.phases()[phaseIndexes_[phasei]].name() == phaseName)
        {
            return true;
        }
    }
    return false;
}


void Foam::kineticTheorySystem::correct()
{
    alphap_ = 0.0;
    Up_ = dimensionedVector("0", dimVelocity, Zero);
    Thetap_ = dimensionedScalar("0", sqr(dimVelocity), 0.0);

    forAll(phaseIndexes_, phasei)
    {
        const phaseModel& phase = fluid_.phases()[phaseIndexes_[phasei]];
        const volScalarField& alpha = phase;

        alphap_ += alpha;
        Up_ += alpha*phase.U();
        Thetap_ +=  alpha*Thetas_[phasei]();
    }
    Up_ /= max(alphap_, residualAlpha_);
    Thetap_ /= max(alphap_, residualAlpha_);

    forAll(pairs_, pairi)
    {
        const phasePairKey& key = pairs_[pairi];
        const phaseModel& phase1 = fluid_.phases()[key.first()];
        const phaseModel& phase2 = fluid_.phases()[key.second()];
        const volScalarField& Theta1 = Thetas_[phase1.name()];
        const volScalarField& Theta2 = Thetas_[phase2.name()];

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

    Pfr_ = frictionalStressModel_->frictionalPressure
    (
        alphap_,
        alphaMax_
    );

    PfrPrime_ = frictionalStressModel_->frictionalPressurePrime
    (
        alphap_,
        alphaMax_
    );

    tmp<volTensorField> tgradU(fvc::grad(Up_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    nuFric_ = dimensionedScalar("0", nuFric_.dimensions(), 0.0);
    forAll(phaseIndexes_, phasei)
    {
        const phaseModel& phase = fluid_.phases()[phaseIndexes_[phasei]];
        nuFric_ += frictionalStressModel_->nu
        (
            phase,
            alphap_,
            alphaMax_,
            Pfr_/phase.rho(),
            D
        )*phase;
    }
    nuFric_ /= max(alphap_, residualAlpha_);

    alphaMax_ = packingLimitModel_->alphaMax();
    alphaMax_.correctBoundaryConditions();
}


void Foam::kineticTheorySystem::correctAlphap()
{
    alphap_ = 0.0;
    forAll(phaseIndexes_, phasei)
    {
        alphap_ += fluid_.phases()[phaseIndexes_[phasei]];
    }
    alphap_.correctBoundaryConditions();
}
// ************************************************************************* //
