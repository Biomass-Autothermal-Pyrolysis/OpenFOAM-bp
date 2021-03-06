/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "PhaseTransferPhaseSystem.H"
#include "phaseTransferModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::addToSpeciesDmdt
(
    const phasePair& pair,
    const wordList species
)
{
    if (species.empty())
    {
        return;
    }

    if (!specieDmdt_.found(pair))
    {
        specieDmdt_.insert(pair, new HashPtrTable<volScalarField>());
    }

    forAllConstIter
    (
        wordList,
        species,
        sIter
    )
    {
        const word specieName(*sIter);
        if (!specieDmdt_[pair]->found(specieName))
        {
            specieDmdt_[pair]->insert
            (
                specieName,
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            specieName+":dmdt",
                            pair.name()
                        ),
                        this->mesh().time().timeName(),
                        this->mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    this->mesh(),
                    dimensionedScalar(Foam::phaseTransferModel::dimDmdt, 0)
                )
            );
        }
    }
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::dmdtTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::
sumDmdts() const
{
    autoPtr<phaseSystem::dmdtTable> sumDmdtsPtr(new phaseSystem::dmdtTable);

    phaseSystem::dmdtTable& sumDmdts = sumDmdtsPtr();

    forAllConstIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[phaseTransferModelIter.key()];

        sumDmdts.insert
        (
            pair,
            zeroVolField<scalar>(pair, "dmdt", dimDensity/dimTime).ptr()
        );

        if (specieDmdt_.found(pair))
        {
            forAllConstIter
            (
                HashPtrTable<volScalarField>,
                *specieDmdt_[pair],
                sDmdtIter
            )
            {
                *sumDmdts[pair] +=
                    phaseTransferModelIter()->speciesDmdt(sDmdtIter.key());
            }
        }
        else
        {
            *sumDmdts[pair] += *rDmdt_[pair];
        }
    }

    return sumDmdtsPtr;
}


// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::addDmdtY
(
    const phaseSystem::dmdtTable& dmdts,
    phaseSystem::specieTransferTable& eqns
) const
{
    forAllConstIter(phaseSystem::dmdtTable, dmdts, dmdtIter)
    {
        const phasePairKey& key = dmdtIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const volScalarField dmdt(Pair<word>::compare(pair, key)**dmdtIter());
        const volScalarField dmdt12(negPart(dmdt));
        const volScalarField dmdt21(posPart(dmdt));

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        // Note that the phase YiEqn does not contain a continuity error term,
        // so the transfers below are complete.

        forAll(phase1.Y(), Yi1)
        {
            const volScalarField& Y1 = phase1.Y()[Yi1];
            const volScalarField& Y2 = phase2.Y(Y1.member());

            *eqns[Y1.name()] += dmdt21*Y2 + fvm::Sp(dmdt12, Y1);
            *eqns[Y2.name()] -= dmdt12*Y1 + fvm::Sp(dmdt21, Y2);
        }
    }
}


template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::addDmidtY
(
    const phaseSystem::dmidtTable& dmidts,
    phaseSystem::specieTransferTable& eqns
) const
{
    forAllConstIter(phaseSystem::dmidtTable, dmidts, dmidtIter)
    {
        const phasePairKey& key = dmidtIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();

        // Note that the phase YiEqn does not contain a continuity error term,
        // so the transfers below are complete.

        forAllConstIter(HashPtrTable<volScalarField>, *dmidtIter(), dmidtJter)
        {
            const word& member = dmidtJter.key();

            const volScalarField dmidt
            (
                Pair<word>::compare(pair, key)**dmidtJter()
            );

            if (!phase1.pure())
            {
                const volScalarField& Y1 = phase1.Y(member);
                *eqns[Y1.name()] += dmidt;
            }

            if (!phase2.pure())
            {
                const volScalarField& Y2 = phase2.Y(member);
                *eqns[Y2.name()] -= dmidt;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::PhaseTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generatePairsAndSubModels
    (
        "phaseTransfer",
        phaseTransferModels_,
        false
    );

    forAllConstIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[phaseTransferModelIter.key()];

        addToSpeciesDmdt(pair, phaseTransferModelIter()->activeSpecies());

        if (!specieDmdt_.found(pair))
        {
            this->rDmdt_.insert
            (
                phaseTransferModelIter.key(),
                zeroVolField<scalar>(pair, "rDmdt", dimDensity/dimTime).ptr()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::
~PhaseTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    autoPtr<phaseSystem::dmdtTable> sumDmdtsPtr = this->sumDmdts();

    const phaseSystem::dmdtTable& sumDmdts = sumDmdtsPtr();

    forAllConstIter(phaseSystem::dmdtTable, sumDmdts, sumDmdtIter)
    {
        const phasePair& pair = this->phasePairs_[sumDmdtIter.key()];

        addField(pair.phase1(), "dmdt", *sumDmdtIter(), dmdts);
        addField(pair.phase2(), "dmdt", - *sumDmdtIter(), dmdts);        

    }

    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtU(sumDmdts(), eqns);
    
    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtU(sumDmdts(), eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::heatTransfer() const
{
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    this->addDmdtHe(rDmdt_, eqns);
    this->addDmidtHe(specieDmdt_, eqns);


    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::specieTransfer() const
{
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr
    (
        new phaseSystem::specieTransferTable()
    );

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    // Create a mass transfer matrix for each species of each phase
    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        const PtrList<volScalarField>& Yi = phase.Y();

        forAll(Yi, i)
        {
            eqns.insert
            (
                Yi[i].name(),
                new fvScalarMatrix(Yi[i], dimMass/dimTime)
            );
        }
    }

    this->addDmdtY(rDmdt_, eqns);
    this->addDmidtY(specieDmdt_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    forAllConstIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
	if (specieDmdt_.found(phaseTransferModelIter.key()))
        {
            forAllIter
            (
                HashPtrTable<volScalarField>,
                *specieDmdt_(phaseTransferModelIter.key()),
                sDmdtIter
            )
            {
                (*sDmdtIter)->ref() = Zero;
            }
        }
        else
        {
            *rDmdt_[phaseTransferModelIter.key()] =
                dimensionedScalar(dimDensity/dimTime, 0);
        }
    }

    forAllIter
    (
        phaseTransferModelTable,
        phaseTransferModels_,
        phaseTransferModelIter
    )
    {
	phaseTransferModelIter()->correct();

        if (specieDmdt_.found(phaseTransferModelIter.key()))
        {
            forAllIter
            (
                HashPtrTable<volScalarField>,
                *specieDmdt_(phaseTransferModelIter.key()),
                sDmdtIter
            )
            {
                (*sDmdtIter)->ref() +=
                    phaseTransferModelIter()->speciesDmdt(sDmdtIter.key());
            }
        }
        else
        {
            *rDmdt_[phaseTransferModelIter.key()] +=
                phaseTransferModelIter()->dmdt();
        }
    }
}


template<class BasePhaseSystem>
bool Foam::PhaseTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
