/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Application
    buoyantBoussinesqPimpleFoam

Description
    Transient solver for turbulent incompressible flow with particle tracking 
    based on icoCoupledKinematicParcelFoam and buoyantBoussinesqPimpleFoam.

    The particle tracking implementation requires a absolute pressure in Pa,
    including the hydrostatic pressure contribution, refered to as 'p'.  The
    turbulence modelling requires a pressure normalised with the densisty
    without the hydrostatic pressure contribution, refered to a 'p_rgh' (just
    as in buoyantBoussinesqPimpleFoam). Therefore a similar split up between
    dynamic and absolute pressure as in buoyantBoussinesqPimpleFoam is made.
    The absolute pressure is multiplied by the density in order to use the
    particle tracking routines as done in icoCoupledKinematicParcelFoam.

    In order to allow for particle tracking with turbulent fluctuations, 
    (set in kinematicCloudProperties with stochasticDispersionRAS), a modified
    version of the lagrangian particle tracking library is used, in which 
    all compressible turbulent models are replaced by incompressible model.


    Eelco van Vliet
    Feb 2012

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "basicKinematicCollidingCloud.H"
#include "pimpleControl.H"
#include "IObasicSourceList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "cloudName",
        "name",
        "specify alternative cloud name. default is 'kinematicCloud'"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "Evolving " << kinematicCloud.name() << endl;
        kinematicCloud.evolve();
        kinematicCloud.info();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
