#/usr/bin/env python


### Units for systems in SIM suite.
### coded by MSS

import numpy as np


class UnitsClass:
    
    ModuleConsts = ["kB", "FPE", "ECharge", "eV"]

    def __init__(self, UnitSystem = {}):
        for (k, v) in UnitSystem.items():
            setattr(self, k, v)

    def __repr__(self):
        return self.Name
        
    def ConvertEne(self, Ene, Temp):
        """Converts energies reported as string with kT into actual numbers; otherwise returns original."""
        if not Ene is None:
            if type(Ene) is str:
                if "kT" in Ene:
                    kT = Temp * self.kB
                if Ene.endswith("kT"):
                    Ene = float(Ene[:-2]) * kT
                elif Ene.endswith("kTperA"):
                    Ene = float(Ene[:-6]) * kT / self.Angstrom
                else:
                    raise ValueError("'%s' not recognized in terms of units.") 
        return Ene
        


MKSUnits = UnitsClass({
            "Name" : "MKSUnits",
            "kB" : 1.380658e-23,        #boltzmann constant in MKS
            "FPE" : 4.334489e27,	    #4piE0 in e^2/m*J
            "ECharge" : 1.60217733e-19, #charge of electron
            "eV" : 1.602176462e-19,		#electron volt
            "Angstrom" : 1.e-10,    	#defines length scale
            "LScale" : 1.e-10,          #characteristic length scale
            "EScale" : 4.e-21,          #characteristic energy scale
            "MScale" : 1.661e-27,       #characteristic mass scale
            "TimeScale" : 1.e-12,       #characteristic time scale -- number of time units in a ps
            "TempScale" : 4.e-21 / 1.380658e-23,  #characteristic temp scale
            "LLabel" : r"$m$",            #length unit label
            "ELabel" : r"$J$",            #energy unit label
            "MLabel" : r"$kg$",           #mass unit label
            "TimeLabel" : r"$s$",         #time unit label
            "TempLabel" : r"$K$",         #temperature unit label
            "TempDefault" : 300.0,
            })

DimensionlessUnits = UnitsClass({
                      "Name" : "DimensionlessUnits",
                      "kB" : 1.,
                      "FPE" : 1.,
                      "ECharge" : 1.,
                      "eV" : 1.,
                      "Angstrom" : 1.,
                      "LScale" : 1.,
                      "EScale" : 1.,
                      "MScale" : 1.,
                      "TimeScale" : 1.,
                      "TempScale" : 1.,
                      "LLabel" : r"$\sigma$",            
                      "ELabel" : r"$\epsilon$",  
                      "MLabel" : r"$m$",   
                      "TimeLabel" : r"$\sigma\sqrt(m/\epsilon)$",  
                      "TempLabel" : r"$\epsilon/k_B$",   
                      "TempDefault" : 1.0,
                      })

AtomicUnits = UnitsClass({
               "Name" : "AtomicUnits",
               "kB" : 0.0019858775,         #boltzmann constant in kcal/mol/K
               "FPE" : 0.0030114705,        #4piE0 in e^2*mol/A*kcal
               "ECharge" : 1.,   			#charge of electron
               "eV" : 23.052347773, 		#electron volt in kcal/mol
               "Angstrom" : 1.,
               "LScale" : 1.,
               "EScale" : 1.,               #in units of kcal/mol
               "MScale" : 1.,
               "TimeScale" : 20.5,          #time scale -- number of time units in a ps
               "TempScale" : 1. / 0.0019858775, #temp scale in K
               "LLabel" : r"$\AA$",            
               "ELabel" : r"$kcal/mol$",  
               "MLabel" : r"$Da$",   
               "TimeLabel" : r"$ps$",  
               "TempLabel" : r"$K$", 
               "TempDefault" : 300.,               
               })    

#fundamental unit definitions for length, energy, mass, charge, and temp
UnitL = np.array([1,0,0,0,0], int)
UnitE = np.array([0,1,0,0,0], int)
UnitM = np.array([0,0,1,0,0], int)
UnitC = np.array([0,0,0,1,0], int)
UnitT = np.array([0,0,0,0,1], int)




   

