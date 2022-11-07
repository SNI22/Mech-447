# =============================================================================
# Thermo Properties
# =============================================================================

# Reading Polynomial parameters and Calculating Thermodynamic parameters

# The thermo class given will be used to calculate the thermodynamic parameters. 
# Also, you'll need [thermoDataNASA-9.yaml]. This file contains two arrays for 
# each species for two separate temperature ranges: a_lo is for  200<𝑇<1000  K 
# and a_hi is for  1000<𝑇<6000  K. We have followed Cantera nomenclature 
# (http://www.cantera.org/docs/sphinx/html/cti/species.html#the-nasa-9-coefficient-polynomial-parameterization).

# 𝐶0𝑝(𝑇)𝑅𝑢=𝑎0𝑇−2+𝑎1𝑇−1+𝑎2+𝑎3𝑇+𝑎4𝑇2+𝑎5𝑇3+𝑎6𝑇4
 
# 𝐻0(𝑇)𝑅𝑢𝑇=−𝑎0𝑇−2+𝑎1𝑙𝑛𝑇𝑇+𝑎2+𝑎32𝑇+𝑎43𝑇2+𝑎54𝑇3+𝑎65𝑇4+𝑎7𝑇
 
# 𝑠0(𝑇)𝑅𝑢=−𝑎02𝑇−2−𝑎1𝑇−1+𝑎2𝑙𝑛𝑇+𝑎3𝑇+𝑎42𝑇2+𝑎53𝑇3+𝑎64𝑇4+𝑎8
 
# where  𝑅𝑢=8.31446𝐽/𝑚𝑜𝑙/𝐾  is the universal gas constant.

# Note that Gordon and McBride start polynomial coefficient numbering at 1 and 
# the last coeffients for enthalpy and entropy are named respectively  𝑎7=𝑏1  
# and  𝑎8=𝑏2 . 
# (http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20020085330_2002141071.pdf)

# =============================================================================
# Notes
# =============================================================================

# Make sure that the file "thermoDataNASA-9.yaml" is on the same root than this 
# script

# Careful, you might need the following function the first time you are using 
# the code: "pip install pyyaml"


import yaml
import numpy as np

class thermo:
    
    def __init__(self, species, MW) :
        """
        species: input string name of species in thermoData.yaml
        M: input (species molecular weight, kg/kmol)
        """
        
        self.Rgas = 8.31446      # J/mol*K
        self.M    = MW
    
        with open("thermoDataNASA-9.yaml") as yfile:
           yfile = yaml.safe_load(yfile)
        
        self.a_lo = yfile[species]["a_lo"]
        self.a_hi = yfile[species]["a_hi"]
        
        self.T_lo = 200.
        self.T_mid = 1000.
        self.T_hi = 6000.
        
    def cymalalo(self) :
        return self.a_lo

    def cymalahi(self) :
        return self.a_hi

    def cp_mole(self,T) :
        """
        return calorific value at cst p in units J/kmol/K
        T: input (K)
        """
        if T<=self.T_mid and T>=self.T_lo :
            a = self.a_lo
        elif T>self.T_mid and T<=self.T_hi :
            a = self.a_hi
        else :
            print ("ERROR: temperature is out of range")

        cp = a[0]/T**2 + a[1]/T + a[2] + a[3]*T + a[4]*T**2.0 + a[5]*T**3 + a[6]*T**4
        
        return cp * self.Rgas
        
    #--------------------------------------------------------

    def cp_mass(self,T) :
        """
        return calorific value at cst P in units of J/kg/K
        T: input (K)
        """
        return self.cp_mole(T)/self.M

    #--------------------------------------------------------
    
    def h_mole(self,T) :
        """
        return enthalpy in units of J/mol
        T: input (K)
        """
        if T<=self.T_mid and T>=self.T_lo :
            a = self.a_lo
        elif T>self.T_mid and T<=self.T_hi :
            a = self.a_hi
        else :
            print ("ERROR: temperature is out of range")

        hrt = -a[0]/T**2 + a[1]*np.log(T)/T + a[2] + a[3]/2*T + a[4]/3*T**2.0 + a[5]/4*T**3 + a[6]/5*T**4 + a[7]/T
        
        return hrt * self.Rgas * T
        
    #--------------------------------------------------------

    def h_mass(self,T) :
        """
        return enthalpy in units of J/kg
        T: input (K)
        """
        return self.h_mole(T)/self.M

    #--------------------------------------------------------
        
    def s_mole(self,T) :
        """
        return entropy in units of J/mol/K
        T: input (K)
        """
        if T<=self.T_mid and T>=self.T_lo :
            a = self.a_lo
        elif T>self.T_mid and T<=self.T_hi :
            a = self.a_hi
        else :
            print ("ERROR: temperature is out of range")
        
        sr = -a[0]/2/T**2 - a[1]/T + a[2]*np.log(T) + a[3]*T + a[4]/2.0*T**2.0 + a[5]/3.0*T**3.0+ a[6]/4.0*T**4.0+ a[8]
        
        return sr * self.Rgas
        
    #--------------------------------------------------------

    def s_mass(self,T) :
        """
        return entropy in units of J/kg/K
        T: input (K)
        """
        return self.s_mole(T)/self.M
    

# =============================================================================
# Test of the code
# =============================================================================

# Return properties of O2
Molecule = "O2"    
MW_molecule = 32               # kg/kmol

t = thermo(Molecule,MW_molecule) # thermo object;
Temperature = 298       #K
print('Entropy of ',Molecule,' at ',str(Temperature), 'K: ',t.s_mole(Temperature),' J/mol')
print('Enthalpy of ',Molecule,' at ',str(Temperature), 'K: ',t.h_mole(Temperature),' J/mol')
print('Specific heat of ',Molecule,' at ',str(Temperature), 'K: ',t.cp_mole(Temperature),' J/mol')

# =============================================================================
# Comments
# =============================================================================

# If you wish to use this piece of code embedded in another code, you might need the following lines:
from ThermoPropertiesNASA import thermo
# Try it out in a new script, if it works, it should display the same results than the previous results
# If it works, just copy the content of the previous section "Test" in your new script, and adapt it to your needs