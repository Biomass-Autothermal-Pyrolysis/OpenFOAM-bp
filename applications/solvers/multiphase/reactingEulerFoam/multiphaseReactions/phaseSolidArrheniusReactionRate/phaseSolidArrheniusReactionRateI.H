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

#include "phaseSolidArrheniusReactionRate.H"
#include "phaseSystem.H"
#include "diameterModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::phaseSolidArrheniusReactionRate::
phaseSolidArrheniusReactionRate
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dictionary& dict
)
:
    ArrheniusReactionRate(species, dict),
    phaseName_(dict.lookup("phase")),
    ob_(ob),
    reactantName_(dict.lookup("reactant")),
    tA_(ob_.lookupObject<volScalarField>(reactantName_)) // maybe define this as a scalar?
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::phaseSolidArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusReactionRate::operator()(p, T, c, li)*tA_()[li];
}


inline Foam::scalar Foam::phaseSolidArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    return ArrheniusReactionRate::ddT(p, T, c, li)*tA_()[li];
}


inline void Foam::phaseSolidArrheniusReactionRate::write(Ostream& os) const
{
    ArrheniusReactionRate::write(os);
    writeEntry(os, "phase", phaseName_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const phaseSolidArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
