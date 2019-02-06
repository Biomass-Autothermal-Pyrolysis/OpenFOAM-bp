/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::kineticTheoryModel

Description
    Kinetic theory particle phase RAS model

    Reference:
    \verbatim
        van Wachem, B. G. M. (2000).
        Derivation, implementation, and validation of computer simulation models
        for gas-solid fluidized beds.
        PhD Thesis, TU Delft.
    \endverbatim

    There are no default model coefficients.

SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "RASModel.H"
#include "eddyViscosity.H"
#include "phaseCompressibleTurbulenceModel.H"
#include "EddyDiffusivity.H"
#include "phaseModel.H"
#include "phasePairKey.H"
#include "dragModel.H"
#include "kineticTheorySystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
:
    public eddyViscosity
    <
        RASModel<EddyDiffusivity<phaseCompressibleTurbulenceModel>>
    >,
    public tmp<kineticTheoryModel>::refCount
{
    // Private data

            typedef
                HashTable
                <
                    autoPtr<phasePair>,
                    phasePairKey,
                    phasePairKey::hash
                >
                phasePairTable;

            //- Const refrence to phase
            const phaseModel& phase_;

            //- Pointer to central kineticTheory model
            kineticTheorySystem& kineticTheorySystem_;


        // Kinetic Theory Model coefficients

            //- Use equilibrium approximation: generation == dissipation
            Switch equilibrium_;

            //- Residual phase fraction
            dimensionedScalar residualAlpha_;

            //- Maximum turbulent viscosity
            dimensionedScalar maxNut_;


        // Kinetic Theory Model Fields

            //- The granular energy/temperature
            volScalarField Theta_;

            //- The granular bulk viscosity
            volScalarField lambda_;

            //- The granular "thermal" conductivity
            volScalarField kappa_;

            //- Frictional viscosity
            volScalarField nuFric_;



    // Private Member Functions

        void correctNut()
        {}

        //- Disallow default bitwise copy construct
        kineticTheoryModel(const kineticTheoryModel&);

        //- Disallow default bitwise assignment
        void operator=(const kineticTheoryModel&);

        //- Lookup or construct kineticTheorySystem
        kineticTheorySystem& lookupOrConstruct(const fvMesh& mesh, const char*) const;


public:

    //- Runtime type information
    TypeName("kineticTheory");


    // Constructors

        //- Construct from components
        kineticTheoryModel
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const phaseModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kineticTheoryModel();


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const
        {
            return this->nut(patchi);
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<volScalarField> pPrime() const;

        //- Return the face-phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<surfaceScalarField> pPrimef() const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Solve the kinetic theory equations and correct the viscosity
        virtual void correct();

        //- Return cont reference to phase model
        const phaseModel& phase() const;

        //- Return const reference to granular temperature field
        const volScalarField& Theta() const;

        //- Return non-const reference to granular temperature field
        volScalarField& Theta();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
