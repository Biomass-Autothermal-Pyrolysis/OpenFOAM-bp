/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
<<<<<<< HEAD:applications/solvers/multiphase/reactingEulerFoam/kineticTheoryModels/granularPressureModel/granularPressureModel/granularPressureModel.H
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
=======
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
>>>>>>> 3431c1d21e2b51a4533243ab4a5f3479adae96b2:applications/solvers/multiphase/reactingEulerFoam/reactingTwoPhaseEulerFoam/twoPhaseCompressibleTurbulenceModels/kineticTheoryModels/granularPressureModel/granularPressureModel/granularPressureModel.H
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
    Foam::kineticTheoryModels::granularPressureModel

SourceFiles
    granularPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef granularPressureModel_H
#define granularPressureModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "kineticTheorySystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class granularPressureModel Declaration
\*---------------------------------------------------------------------------*/

class granularPressureModel
{
protected:

    // Protected data

        const dictionary& dict_;

        const kineticTheorySystem& kt_;


public:

    //- Runtime type information
    TypeName("granularPressureModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        granularPressureModel,
        dictionary,
        (
            const dictionary& dict,
            const kineticTheorySystem& kt
        ),
        (dict, kt)
    );


    // Constructors

        //- Construct from components
        granularPressureModel
        (
            const dictionary& dict,
            const kineticTheorySystem& kt
        );

        //- Disallow default bitwise copy construction
        granularPressureModel(const granularPressureModel&) = delete;


    // Selectors

        static autoPtr<granularPressureModel> New
        (
            const dictionary& dict,
            const kineticTheorySystem& kt
        );


    //- Destructor
    virtual ~granularPressureModel();


    // Member Functions

        //- Granular pressure coefficient
        virtual tmp<volScalarField> granularPressureCoeff
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& Theta1,
            const volScalarField& Theta2,
            const volScalarField& g0,
            const dimensionedScalar& e
        ) const = 0;

        //- Derivative of the granular pressure coefficient
        virtual tmp<volScalarField> granularPressureCoeffPrime
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const volScalarField& Theta1,
            const volScalarField& Theta2,
            const volScalarField& g0,
            const volScalarField& g0prime,
            const dimensionedScalar& e
        ) const = 0;

        virtual bool read()
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const granularPressureModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
