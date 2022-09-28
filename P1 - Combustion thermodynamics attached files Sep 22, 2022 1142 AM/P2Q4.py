from ipynb.fs.full.ThermoPropertiesNASA import thermo
from P2Q3 import phi #h_gap_rich, R #h_gap_lean, R
from scipy.integrate import quad 
from scipy.optimize import fsolve as solve
import numpy as np

Tf_int = np.empty(21)
T_ref = int(300)

def Cp_lean(Tad):
    return 8 * phi[a] * thermo('CO2',44).cp_mole(Tad) + 9 * phi[a] * thermo('H2O',18).cp_mole(Tad) \
        + (12.5 - 12.5 * phi[a]) * thermo('O2',32).cp_mole(Tad) + 47 * thermo('N2',28).cp_mole(Tad)

def h_int_lean(Tad):
    #print(quad(Cp_lean(Tad), T_ref, Tad)[0])
    return quad(Cp_lean, T_ref, Tad)[0]

def h_gap_lean():
    x = 0.0
    x = phi[a]
    delta_h = - x * thermo('C8H18',114).h_mole(T_ref) - 12.5 * thermo('O2',32).h_mole(T_ref) \
        + 8 * x * thermo('CO2',44).h_mole(T_ref) + 9 * x * thermo('H2O',18).h_mole(T_ref) \
        + (12.5 - 12.5 * x) * thermo('O2',32).h_mole(T_ref) 
    return delta_h

def lean_gap(Tad):
    return h_gap_lean() + h_int_lean(Tad)


def Cp_rich(Tad):
    Cp = 0.8 * phi[a] * thermo('CO',28).cp_mole(Tad) \
    + 7.2 * phi[a] * thermo('CO2',44).cp_mole(Tad) \
    + (25-15.2*phi[a]) * thermo('H2O',18).cp_mole(Tad) \
    + (24.2 * phi[a] - 25) * thermo('H2',2).cp_mole(Tad) + 47 * thermo('N2',28).cp_mole(Tad)

    return Cp

def h_int_rich(Tad):
    return quad(Cp_rich, T_ref, Tad)[0]

def h_gap_rich():
    delta_h = - phi[a] * thermo('C8H18',114).h_mole(T_ref) - 12.5 * thermo('O2',32).h_mole(T_ref) \
    + 0.8 * phi[a] * thermo('CO',28).h_mole(T_ref) \
    + 7.2 * phi[a] * thermo('CO2',44).h_mole(T_ref) \
    + (25-15.2*phi[a]) * thermo('H2O',18).h_mole(T_ref) \
    + (24.2 * phi[a] - 25) * thermo('H2',2).h_mole(T_ref) 
    return delta_h

def rich_gap(Tad):
    return h_gap_rich() + h_int_rich(Tad)

    

for a in range(len(phi)):
    if phi[a] == 1 or phi[a] < 1:
        Tf_int[a] = solve(lean_gap, 1000)
        print(Tf_int[a])
    else:
        Tf_int[a] = solve(rich_gap, 1000)
        print(Tf_int[a])
print(Tf_int)