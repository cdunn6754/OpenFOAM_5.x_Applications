// Here we want to calculate the soot volume fraction as is done in the 
// radiation models.


basicSpecieMixture& mixture_ = thermo.composition();

const volScalarField& T = thermo.T(); 
const volScalarField& p = thermo.p();

scalarField sootSpecificMixtureVolume(T.primitiveField().size(), 0.0);
scalarField rhoInv(T.primitiveField().size(),0.0);

// Sum over all specific mixture volumes 
// units of that would be [Vol_i / mass_mixture]
forAll(mixture_.Y(), specieI)
{
    const scalarField& Yi = mixture_.Y()[specieI];
    const word specieName = mixture_.Y()[specieI].name();
    
    // If it's soot dont use IGL. Use known and constant soot density
    if (specieName == "SOOT")
    {
        forAll(rhoInv, iCell)
        {
            rhoInv[iCell] +=
                Yi[iCell]/rho_soot;
        }
    }
    else // Just use the IGL to find density.
    {
        forAll(rhoInv, iCell)
        {
            rhoInv[iCell] +=
                Yi[iCell]/mixture_.rho(specieI, p[iCell], T[iCell]);
        }   
    }
}
const scalarField& Yj = mixture_.Y("SOOT");

// Now find the soot one once again
sootSpecificMixtureVolume = Yj/rho_soot;

// divide by the summed values to find volume fraction 
// [Vol_i / Vol_mixture]
Vsoot.primitiveFieldRef() = sootSpecificMixtureVolume/rhoInv;
