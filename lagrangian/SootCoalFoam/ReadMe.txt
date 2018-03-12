12-19-17:

This is a solver based very closely on coalChemistryFoam. It uses the PcReactingMultiphase library which enables the secondary pyrolysis rates from PC Coal Lab to be used directly. It also includes (hard coded into the solver) an euler implicit solution scheme for the two equation soot model (soot mass fraction and particle number density) taken from Kronenburg CMC paper. 
