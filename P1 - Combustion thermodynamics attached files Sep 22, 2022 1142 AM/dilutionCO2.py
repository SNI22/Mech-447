"""
Adiabatic flame temperature and equilibrium composition for a fuel/air mixture
as a function of equivalence ratio, including formation of solid carbon.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
Keywords: equilibrium, combustion, multiphase
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import sys
import csv

##############################################################################
# Edit these parameters to change the initial temperature, the pressure, and
# the phases in the mixture.

T = 800
P = 30*101325.0
phi = 0.9

# phases
gas = ct.Solution('Jerzembeck.cti')
carbon = ct.Solution('Jerzembeck.cti')

# the phases that will be included in the calculation, and their initial moles
mix_phases = [(gas, 1.0), (carbon, 0.0)]

# gaseous fuel species
fuel_species = 'IXC8H18'

# equivalence ratio range
npoints = 20

##############################################################################

mix = ct.Mixture(mix_phases)

# create some arrays to hold the data
tad = np.zeros(npoints)
xeq = np.zeros((mix.n_species, npoints))

for i in range(npoints):
    # set the gas state
    gas.set_equivalence_ratio(phi, fuel_species, oxidizer='O2:1.0, N2:3.76', diluent="CO2", fraction={"diluent":(0.05*i)})

    # create a mixture of 1 mole of gas, and 0 moles of solid carbon.
    mix = ct.Mixture(mix_phases)
    mix.T = T
    mix.P = P

    # equilibrate the mixture adiabatically at constant P
    mix.equilibrate('HP', solver='auto', max_steps=1000)

    tad[i] = mix.T
    print('At mole fraction = {0:12.4g}, Tad = {1:12.4g}'.format(0.05*i, tad[i]))
    xeq[:, i] = mix.species_moles

#comp_O = xeq[2,:]
comp_CO = xeq[12,:]
comp_CO2 = xeq[14,:]
#comp_H2 = xeq[4,:]
#comp_H2O = xeq[8,:]
#comp_H = xeq[5,:]
#comp_OH = xeq[6,:]
comp_NO = xeq[104,:]
comp_NO2 = xeq[105,:]

# Define Data

x = np.zeros(npoints)
for i in range(npoints):
    x[i] = 0.05*i
  
# Create Plot

fig, ax1 = plt.subplots() 
  
ax1.set_xlabel('Mole Fraction of CO2') 
ax1.set_ylabel('T_ad', color = 'black') 
ax1.plot(x, tad, color = 'black') 
ax1.tick_params(axis ='y', labelcolor = 'black') 
  
# Adding Twin Axes

ax2 = ax1.twinx() 

ax2.set_yscale('log')
ax2.set_ylabel('Equilibrium Composition', color = 'blue') 
ax2.plot(x, comp_NO, label = 'NO') 
ax2.plot(x, comp_CO, label = 'CO') 
ax2.plot(x, comp_CO2, label = 'CO2') 
ax2.plot(x, comp_NO2, label = 'NO2') 
ax2.tick_params(axis ='y', labelcolor = 'blue') 
 
# Show plot
plt.legend()
plt.show()

# write output CSV file for importing into Excel
csv_file = 'dilution_CO2.csv'
with open(csv_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Mole Fraction of CO2', 'T (K)'] + mix.species_names)
    for i in range(npoints):
        writer.writerow([0.05*i, tad[i]] + list(xeq[:, i]))
print('Output written to {0}'.format(csv_file))

#if '--plot' in sys.argv:
#    import matplotlib.pyplot as plt
#    plt.plot(phi, tad)
#    plt.xlabel('Equivalence ratio')
#    plt.ylabel('Adiabatic flame temperature [K]')
#    plt.show()