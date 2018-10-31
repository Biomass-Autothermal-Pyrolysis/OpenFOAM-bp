/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fv::meanPhaseVelocityForce

Description
    Calculates and applies the force necessary to maintain the specified mean
    velocity.

    Note: Currently only handles kinematic pressure (incompressible solvers).

Usage
    Example usage:
    \verbatim
    selectionMode   all;                    // Apply force to all cells
    fields          (U);                    // Name of velocity field
    phaseName       (air)                   // Name of phase
    Ubar            (10.0 0 0);             // Desired mean velocity
    relaxation      0.2;                    // Optional relaxation factor
    \endverbatim

SourceFiles
    meanPhaseVelocityForce.C

\*---------------------------------------------------------------------------*/

#ifndef meanPhaseVelocityForce_H
#define meanPhaseVelocityForce_H

#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class meanPhaseVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class meanPhaseVelocityForce
:
    public cellSetOption
{
protected:

    // Protected data

        //- Phase name
        const word phaseName_;

        //- Constant reference to alpha field
        const volScalarField& alpha_;

        //- Constant reference to density field
        const volScalarField& rho_;

        //- Average velocity
        vector Ubar_;

        //- Pressure gradient before correction
        scalar gradP0_;

        //- Change in pressure gradient
        scalar dGradP_;

        //- Flow direction
        vector flowDir_;

        //- Relaxation factor
        scalar relaxation_;

        //- Matrix 1/A coefficients field pointer
        autoPtr<volScalarField> rAPtr_;


    // Protected Member Functions

        //- Calculate and return the magnitude of the mean velocity
        //  averaged over the selected cellSet
        virtual scalar magUbarAve(const volVectorField& U) const;

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const scalar gradP) const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        meanPhaseVelocityForce(const meanPhaseVelocityForce&);

        //- Disallow default bitwise assignment
        void operator=(const meanPhaseVelocityForce&);


public:

    //- Runtime type information
    TypeName("meanPhaseVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        meanPhaseVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluate

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
