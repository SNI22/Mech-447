from ipynb.fs.full.ThermoPropertiesNASA import thermo
import numpy as np
import matplotlib.pyplot as plt

phi=np.arange(0.5,1.55,0.05)
print(phi[20])
Tf = np.empty(21)

T_ref = 300 #K
R = 8.3144621
Cp = 3.5 * R

def h_gap_lean(phi):
    delta_h = - phi * thermo('C8H18',114).h_mole(T_ref) - 12.5 * thermo('O2',32).h_mole(T_ref) \
    + 8 * phi * thermo('CO2',44).h_mole(T_ref) \
    + 9 * phi * thermo('H2O',18).h_mole(T_ref) \
    + (12.5 - 12.5 * phi) * thermo('O2',32).h_mole(T_ref)
    
    return delta_h

def h_gap_rich(phi):
    delta_h = - phi * thermo('C8H18',114).h_mole(T_ref) - 12.5 * thermo('O2',32).h_mole(T_ref) \
    + 0.8 * phi * thermo('CO',28).h_mole(T_ref) \
    + 7.2 * phi * thermo('CO2',44).h_mole(T_ref) \
    + (25-15.2*phi) * thermo('H2O',18).h_mole(T_ref) \
    + (24.2 * phi - 25) * thermo('H2',2).h_mole(T_ref) 
    
    return delta_h



for a in range (len(phi)):
    if phi[a] == 1 or phi[a] < 1:
        Tf[a] = (- h_gap_lean(phi[a])) /(Cp*(4.5*phi[a] + 47 + 12.5)) + 300
    else:
        Tf[a] = (- h_gap_rich(phi[a])) / (Cp*( - 15.2 * phi[a] + 8*phi[a] + 24.2*phi[a] + 47)) + 300

print(Tf)


