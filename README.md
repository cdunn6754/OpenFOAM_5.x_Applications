# OpenFOAM 5.x Solvers and Applications

## Description of Contents

These solvers have been developed in a large part to include a 
more accurate treatment of soot and coal volatiles. They 
are all built upon the OpenFOAM (OF) solver CoalChemistryFoam. There
are two main branches, those that rely on P.C. Coal Lab (PCCL) for devolatilization
and tar evolution predictions, and those that rely on empirical findings (e.g. from
Ma 1996 and Brown 1997).

### Solvers (within lagrangian):

This folder is called lagrangian because I was attempting to mirror the OF
default file system and all of the solvers we create involve lagrangian modeling
of coal particles.

  * SootCoalFoam
    This solver uses PCCL to predict devolatilization. It utilizes the 
    [PCReactingMultiphase Library](https://github.com/cdunn6754/OpenFOAM_5.x_Libraries/tree/master/lagrangian/PcReactingMultiphase)
    which pretty drastically modifies the underlying cloud and parcel classes. 
    
  * SootTarFoam
    This solver uses empirical rates to describe the decomposition of tar. It uses the TwoEquationSoot 
    and SootTarModel.
    
  * NoGrowthSootTarFoam
    The same as the SootTarFoam solver but uses the NoGrowthTwoEquationSoot model rather than TwoEquationSoot.
    
  * PCFSootTarFoam
    Same as NoGrowthSootTarFoam but uses the PCFSootTarModel and so it 
	  incorporates the [PCF](https://github.com/cdunn6754/PCCLConversion) predicted rate constant.
    
### Utilities:

  * PCCL2Foam
    A handy python script that can be used to convert PCCL outputs into OF inputs for devolatilization predictions.
    
   * sootVolumeFraction
    This is a solver that I have pared down to just be run in a case and calculate the soot volume fraction 
    field. It automatically does this for the latest timestep and writes the new field within the timestep folder.
    
 ## Building:
 
 To build these solvers you will need to have OpenFOAM 5.x installed and have already built the corresponding
 [libraries](https://github.com/cdunn6754/OpenFOAM_5.x_Libraries)  Then you should clone this repo to your
 user application directory, e.g. `~/OpenFOAM/<username>-5.x/ applications/`. Beyond that you can use the OF 
 wmake command to compile the solvers. Each of them has its own Make folder with contents `files` and `options`. 
 If you have an error in compilation the first place to start is always by making sure the links to libraries
 and included header files are accurate within `Make/options`. 
    
