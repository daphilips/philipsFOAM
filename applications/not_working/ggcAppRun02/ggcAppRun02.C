/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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
    ggcAppRun02

Description
    Steady-state solver for incompressible, turbulent flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    // BEGIN set boundary conditions

    // Get index of patch 
    label inletPatchID = mesh.boundaryMesh().findPatchID("inlet_1"); 

    // Get reference to boundary value 
    const fvPatchVectorField& faceCentreshub = mesh.C().boundaryField()[inletPatchID]; 
    fvPatchVectorField& inlet_1U = U.boundaryField()[inletPatchID]; 
    fvPatchScalarField& inlet_1k = k.boundaryField()[inletPatchID];  
 
    //constants for profile calculation
    double u_max = 8.0;
    double u_min = 4.5;
    double y_top = 1.4;
    double y_mid = 0.6;
    
    double u_temp;
    // loop over all hub faces 
    forAll(inlet_1U, faceI) 
    { 
      // get coordinate for face centre 
      const vector& c = faceCentreshub[faceI]; 
      u_temp = (u_max-u_min)*( (c[1]-y_mid)/(y_top-y_mid) ) + u_min;
      if (c[1] > 1.38)
        u_temp = ( 7.9125*(y_top - c[1])/(y_top-1.38) );
             
      vector u_prof(u_temp, 0.0, 0.0);
  
      inlet_1U[faceI] = u_prof; 
    }

    double TI = 0.05;

    forAll(inlet_1k,faceI)
    {
     // get coordinate for face centre
      const vector& c = faceCentreshub[faceI];
      u_temp = (u_max-u_min)*( (c[1]-y_mid)/(y_top-y_mid) ) + u_min;
      
      scalar k_prof(((u_temp*TI)*(u_temp*TI))*3.0/2.0);
    
      inlet_1k[faceI] = k_prof;

    }

    // Force the write 
    U.write(); 
//    k.write();
    
    // END set boundary conditions
    Info<< "\n ExecutionTime = " 
        << runTime.elapsedCpuTime() 
        << " s\n" << endl; 

    Info<< "End" << endl; 

    return 0;
}


// ************************************************************************* //
