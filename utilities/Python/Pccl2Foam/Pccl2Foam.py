# Convert some of the PCCL outputs to usable inputs in OpenFOAM cloud properties files.

### Inputs from proximate and ultimate analysis (not PCCL)
Y_ash = .0845
Y_liq = 0.0758
# assume the liquid is water, then find 
# fraction that is daf
Y_daf = 1 - Y_ash - Y_liq

### Inputs from PCCL all are given on daf basis
# Relevant File name: "FDC1WT1.RPT"
Y_char_daf = 0.402
Y_gas_daf = 0.598 #NOTE: sum both Ygas and Ytar from file

## dict of primary volatile species daf mass fraction
# Only include species which are included in the combustion mechanism
# Don't worry if they dont sum to one we will renormalize later
# Relevant File names:
#   "FDC1WT1.RPT","FDC1HC1.RPT", "FDC1NG1.RPT"
Y_i_daf = {
    "TAR":0.417,
    "CO2":0.027,
    "H2O":0.055,
    "CO":0.027,
    "H2":0.0091,
    "CH4":0.028,
    "C2H2":0.0105,
    "C2H6":0.0103
}

## Secondary Tar Breakdown
# Relevant File Name: "SFTRC1T1.RPT"
sec_species = {
    "SOOT":.602,
    "C2H2":.419e-1,
    "H2":.524e-1,
    "CO":.969e-1,
    "H2O":0.137e-2,
    "CO2":0.119e-1
}

### **** END INPUTS **** ###


### Convert Fractions to a non-daf basis
## Phase Fractions
Y_char = Y_char_daf * Y_daf
Y_solid = Y_char + Y_ash
Y_gas = 1.0 - Y_solid - Y_liq # includes tar too

## Solid Phase Fractions 
## (i.e. fraction of specie within the solid phase as in OpenFOAM dict)
Y_char_solid = Y_char/Y_solid
Y_ash_solid = Y_ash/Y_solid

## Liquid Phase Fractions (for this all liquid is water)
Y_water_liq = 1.0

## Gas Phase Fractions 
# here we make the assumption that Y_i_daf/Y_gas_daf = Y_i/Y_gas
# which I think is pretty much true.
# Note that tar fraction is preserved and the rest are normalized to 
# account for missing species.
gas_sum_unnormalized = sum([Y_i_daf[name] for name in Y_i_daf.keys() if name != "TAR"])
Y_gas_no_tar = Y_gas - Y_i_daf["TAR"]
for species in Y_i_daf.keys():
    if species == "TAR":
        continue
    old_mf = Y_i_daf[species]
    Y_i_daf[species] = (old_mf/gas_sum_unnormalized) * Y_gas_no_tar

# Now we need to normalize the entire array to total 1.0
# i.e. currently sum(Y_i_daf.values()) == Y_gas, but we want
# sum(Y_i_daf.values()) == 1.0
for species in Y_i_daf:
    Y_i_daf[species] =Y_i_daf[species] / Y_gas
    
### Report the phase fractions (gas, liquid, solid) (in OpenFOAM ready form)
print("\nPhase Fractions:\n Gas: {:.3f}\n Liquid: {:.3f}\n Solid: {:.3f}\n".format(
    Y_gas, Y_liq, Y_solid))


### Report the fractions within each phase (in OpenFOAM ready form)
print("\nGas Fractions:")
for specieName in Y_i_daf.keys():
    Yspecie = Y_i_daf[specieName]
    print(" {}: {:.3f}".format(
        specieName, Yspecie))

print("\nLiquid Fractions")
print(" H2O: {:.3f}".format(Y_water_liq))

print("\nSolid Fractions")
print(" Char: {:.3f}".format(Y_char_solid))
print(" Ash: {:.3f}".format(Y_ash_solid))
      

### Normalize and Report Secondary Tar products ratios (in OpenFOAM ready form)
print("\nSecondary Devolatilization:")
for specieName in sec_species.keys():
    Yspecie = sec_species[specieName]
    
    sec_species[specieName] = Yspecie/sum(sec_species.values())
    print(" {}: {:.3f}".format(
        specieName, sec_species[specieName]))



