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

Class
    Foam::phaseSolidArrheniusReactionRate

Description
    A modified Arrhenius reaction rate given by:

        k = (A * T^beta * exp(-Ta/T))*c

    Where c is the concentration, the name of which is
    specified by the user.

SourceFiles
    phaseSolidArrheniusReactionRateI.H
    
\*---------------------------------------------------------------------------*/

#ifndef phaseSolidArrheniusReactionRate_H
#define phaseSolidArrheniusReactionRate_H

#include "ArrheniusReactionRate.H"
#include "speciesTable.H"
#include "phaseSystem.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class phaseSolidArrheniusReactionRate;

Ostream& operator<<(Ostream&, const phaseSolidArrheniusReactionRate&);

/*---------------------------------------------------------------------------*\
             Class surfaceArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class phaseSolidArrheniusReactionRate
:
    public ArrheniusReactionRate
{
    // Private Data

        //- Name of the phase
        const word phaseName_;

        //- Reference to the object registry
        const objectRegistry& ob_;

        //- Reactant name
        const word reactantName_;

        //- concentration of speice
        mutable tmp<volScalarField> tA_;


public:

    // Constructors

        //- Construct from dictionary
        inline phaseSolidArrheniusReactionRate
        (
            const speciesTable& species,
            const objectRegistry& ob,
            const dictionary& dict
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "phaseSolidArrhenius";
        }

        //- Evaluate the rate
        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Evaluate the derivative
        inline scalar ddT
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) const;

        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const phaseSolidArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseSolidArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
