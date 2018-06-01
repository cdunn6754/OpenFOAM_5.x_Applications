'''
Convert some of the PCCL outputs to usable inputs in OpenFOAM cloud properties files.

06-01-18
Python 3.6.1
'''

### User Inputs

# Do we need to analyze the secondary breakdown products too?
analyze_secondary_breakdown = False

### Inputs from proximate and ultimate analysis (not PCCL)
Y_ash = 0.0617
Y_liquid = 0.268
# assume the liquid is water, then find 
# fraction that is daf
Y_daf = 1.0 - Y_ash - Y_liquid

### Inputs from PCCL all are given on daf basis
# Relevant File name: "FDC1WT1.RPT"
Y_char_daf = 0.397
Y_gas_daf = 1.0 - Y_char_daf #NOTE: sum of both Ygas and Ytar from file

## dict of primary volatile species daf mass fraction
# Only include species which are included in the combustion mechanism
# Don't worry if they dont sum to one we will renormalize later
# In this example I have lumped all C2 and C3 HC in the C2 species we
# can handle, that is recommended.
# Relevant File names:
#   "FDC1WT1.RPT","FDC1HC1.RPT", "FDC1NG1.RPT"
Y_i_daf = {
    "TAR":0.298,
    "CO2":0.074,
    "H2O":0.077,
    "CO":0.069,
    "H2":0.0066,
    "CH4":0.041,
    "C2H4":0.0164,
    "C2H6":0.0153   ##Added C3H6 and C2H6 together
}

if analyze_secondary_breakdown:
    ## Secondary Tar Breakdown
    # Relevant File Name: "SFTRC1T1.RPT"
    sec_species = {
        "SOOT":.748, ## If necessary be sure to add Oil and PAH fraction here
        "C2H2":.322e-1,
        "H2":.470e-1,
        "CO":.146,
        "H2O":.566e-3,
        "CO2":.365e-2
    }

### **** END INPUTS **** ###




### Convert Fractions to a non-daf basis
## Phase Fractions
Y_gas = Y_gas_daf * Y_daf # includes tar too
Y_char = 1.0 - Y_gas - Y_liquid - Y_ash
Y_solid = Y_char + Y_ash


## Solid Phase Fractions (fraction of species within phase)
## (i.e. fraction of specie within the solid phase as in OpenFOAM dict)
Y_char_solid = Y_char/Y_solid
Y_ash_solid = Y_ash/Y_solid

## Liquid Phase Fractions (fraction of species within phase)
Y_water_liq = 1.0

## Gas Phase Fractions (fraction of specie within phase)
# here we make the assumption that Y_i_daf/Y_gas_daf = Y_i/Y_gas
# which I think is pretty much true.
# Note that tar fraction is preserved and the rest are normalized to 
# account for missing species.

gas_sum_no_tar = sum([Y_i_daf[name] for name in Y_i_daf.keys() if name != "TAR"])

for species in Y_i_daf.keys():
    if species == "TAR":
        continue
    old_mf = Y_i_daf[species]
    Y_i_daf[species] = ( (old_mf/gas_sum_no_tar)  
                         * (Y_gas_daf - Y_i_daf["TAR"]))

# Species fractions of gas, not daf
Y_i_gas = {}

# Now we need to normalize the entire array to total 1.0
# i.e. currently sum(Y_i_daf.values()) == Y_gas_daf, but we want
# sum(Y_i.values()) == 1.0
for species in Y_i_daf:
    Y_i_gas[species] =Y_i_daf[species] / Y_gas_daf

### END Analysis
    


### Report the mass fraction that is DAF
print("\nDAF mass fraction: \n Y_daf: {}".format(Y_daf))

### Report the phase fractions (gas, liquid, solid) (in OpenFOAM ready form)
print("\nPhase Fractions:\n Gas: {:.3f}\n Liquid: {:.3f}\n Solid: {:.3f}".format(
    Y_gas, Y_liquid, Y_solid))

### Report the fractions within each phase (in OpenFOAM ready form)
print("\nGas Fractions:")
for specieName in Y_i_daf.keys():
    Yspecie = Y_i_gas[specieName]
    print(" {}: {:.3f}".format(
        specieName, Yspecie))

print("\nLiquid Fractions")
print(" H2O: {:.3f}".format(Y_water_liq))

print("\nSolid Fractions")
print(" Char: {:.3f}".format(Y_char_solid))
print(" Ash: {:.3f}".format(Y_ash_solid))
      

if analyze_secondary_breakdown:
    ### Normalize and Report Secondary Tar products ratios (in OpenFOAM ready form)
    print("\nSecondary Devolatilization:")
    for specieName in sec_species.keys():
        Yspecie = sec_species[specieName]
    
        sec_species[specieName] = Yspecie/sum(sec_species.values())
        print(" {}: {:.3f}".format(
            specieName, sec_species[specieName]))

print("\n\n")



