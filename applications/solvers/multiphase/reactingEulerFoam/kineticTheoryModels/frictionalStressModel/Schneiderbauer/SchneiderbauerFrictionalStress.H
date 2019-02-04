/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::kineticTheoryModels::frictionalStressModels::Schneiderbauer

Description

SourceFiles
    SchneiderbauerFrictionalStress.C

\*---------------------------------------------------------------------------*/

#ifndef SchneiderbauerFrictionalStress_H
#define SchneiderbauerFrictionalStress_H

#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{

/*---------------------------------------------------------------------------*\
                           Class Schneiderbauer Declaration
\*---------------------------------------------------------------------------*/

class Schneiderbauer
:
    public frictionalStressModel
{
    // Private data

        dictionary coeffDict_;

        //- Coefficients
        dimensionedScalar phi1_;
        dimensionedScalar phi2_;
        dimensionedScalar Io_

        //- Minimum volume fraction that frictional pressures act
        dimensionedScalar alphaMinFriction_;


public:

        //- Runtime type information
        TypeName("Schneiderbauer");


    // Constructors

        //- Construct from components
        Schneiderbauer(const dictionary& dict);


    //- Destructor
    virtual ~Schneiderbauer();


    // Member functions

        virtual tmp<volScalarField> frictionalPressure
        (
            const volScalarField& alphap,
            const volScalarField& alphaMax
        ) const;

        virtual tmp<volScalarField> frictionalPressurePrime
        (
            const volScalarField& alphap,
            const volScalarField& alphaMax
        ) const;

        virtual tmp<volScalarField> nu
        (
            const phaseModel& phase,
            const volScalarField& alphap,
            const volScalarField& alphaMax,
            const volScalarField& pf,
            const volSymmTensorField& D
        ) const;

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frictionalStressModels
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
