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

Utility
    sootVolumeFraction

Description
    Based on the soot mass fraction and soot particle number density
    determine the particle diameter and the soot volume fraction on a
    cell basis.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// for thermo
#include "psiCombustionModel.H"
#include "turbulentFluidThermoModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

    // Using the equation from Kronenburg CMC paper published 2000
    // for soot surface area and particle diameters
    // density from Dasgupta thesis 2015
    const scalar rho_soot = 2000.0; //[kg/m^3]
    const scalar factor = 6.0/(Foam::constant::mathematical::pi * rho_soot);


    // Use the other method to get a value for sootVolumeFraction
        #include "alternativeVsoot.H"
    
    // loop through cells
    forAll(Dp, cellNumber)
    {
        scalar& cellDp = Dp[cellNumber];
        scalar& cellSootVFraction = Vsoot[cellNumber];
        const scalar& cellNs = Ns[cellNumber];
        const scalar& cellYsoot = Ysoot[cellNumber];
        //const scalar& cellVolume = mesh.V()[cellNumber];
        //const scalar& cellRho = rho[cellNumber];


        // If Ns is too small just forget about this cell
        if (cellNs <= 1e-3)
        {
            cellDp = 0.0;
            cellSootVFraction = 0.0;
        }
        else
        {
            // calculate soot diameters in this cell (uniform diameters per cell) 
            // From Kronenburg equation
            cellDp = Foam::pow(factor * (cellYsoot/cellNs), 1.0/3.0);


            //****************************//
            //    DEPRECATED: this method relies on the density calculation
            //    from thermo that accounts for soot with the IGL. We need
            //    in some cases to use the volume fraction densitiy calulcation.
            //    That calculation is carried out in "alternativeVsoot.H" now.
            //*****************************//

            // // Volume of soot particles in cell 
            // scalar sootParticleV = (1.0/6.0) * Foam::constant::mathematical::pi
            //     * Foam::pow(cellDp,3.0);
            
            // // Ns is the #particles/kg we need to find the mass in the cell
            // // then find total soot particle volume in cell
            // // [=] ([particles/kg] * [kg/m^3] * [m^3] * [m^3(soot)/particle] 
            // scalar cellSootVolume = (cellNs*cellRho)* cellVolume * sootParticleV;
            // cellSootVFraction = cellSootVolume/cellVolume;
        }
    }
    // Write the two new fields to the time folder.
    Info << "Writing Vsoot and Dp fields to time folder " << runTime.timeName()
        << ".\n" << endl;
    Vsoot.write();
    Dp.write();

    Info << "The soot particle diameters are \nmin: " << min(Dp.primitiveField())
        << "\nmax: " << max(Dp.primitiveField()) << "\n" << endl;

    Info << "The soot volume fractions are \nmin: " << min(Vsoot.primitiveField())
        << "\nmax: " << max(Vsoot.primitiveField()) << endl;


    // // Debugging to find the difference in the two methods
    // scalar maxDiff(0.0);
    // label maxDiffNumber(0);
    // scalar minDiff(10000.0);
    // label minDiffNumber(0);
    // scalar cellDiff(0.0);

    // forAll(sootVolumeRad, cellNumber)
    // {
    //     if (Vsoot[cellNumber] > 0.0)
    //     {
    //         cellDiff = (sootVolumeRad[cellNumber]/ Vsoot[cellNumber]);

    //         if (cellDiff > maxDiff)
    //         {
    //             maxDiff = cellDiff;
    //             maxDiffNumber = cellNumber;

    //         }

    //         if (cellDiff < minDiff)
    //         {
    //             minDiff = cellDiff;
    //             minDiffNumber = cellNumber; 
    //         }
    //     }
    // }

    // Info << "\nNew Max Diff: " <<  maxDiff << endl;
    // Info << "Cell Number: " << maxDiffNumber << endl;
    // Info << "RAD: " << sootVolumeRad[maxDiffNumber] << endl;
    // Info << "Original: " << Vsoot[maxDiffNumber] << endl;
    // Info << "Ysoot: " << Ysoot[maxDiffNumber] << endl;
    // Info << "Ns: " << Ns[maxDiffNumber] << endl;
    // Info << "rho: " << rho[maxDiffNumber] << endl;
    // Info << "Volume: " << mesh.V()[maxDiffNumber] << endl;
    

    // Info << "\nNew Min Diff: " <<  minDiff << endl;
    // Info << "Cell Number: " << minDiffNumber << endl;
    // Info << "RAD: " << sootVolumeRad[minDiffNumber] << endl;
    // Info << "Original: " << Vsoot[minDiffNumber] << endl;
    // Info << "Ysoot: " << Ysoot[minDiffNumber] << endl;
    // Info << "Ns: " << Ns[minDiffNumber] << endl;
    // Info << "rho: " << rho[minDiffNumber] << endl;
    // Info << "Volume: " << mesh.V()[minDiffNumber] << endl;

    Info << "\nThe one: " << Vsoot[826];

    return 0;
}


// ************************************************************************* //
