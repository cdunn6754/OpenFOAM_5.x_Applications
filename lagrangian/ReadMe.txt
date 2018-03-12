SootCoalFoam: 
01-30-18:
Solver based on coalChemistryFoam. Built for use with the PcReactingMultiphase library (i.e. can utilize devol rates and secondary tar breakdown rates from PC coal lab). Also built for use with the EulerImplicitSystem library which adds the capability to solve the 2-equation soot models.

sootCoalChemistryFoam:
01-30-18:
Only created to help debug the sootReactingThermo library. We are trying to force the thermo classes to ignore soot in the calculation of gas density. Then we can use the method from the Zimmer (2016) paper to find the total density 
     rho = (f_v)*rho_soot + (1.0 - f_v)*rho_gas.

03-06-18:

Deleted the sootCoalChemistyFoam solver.
