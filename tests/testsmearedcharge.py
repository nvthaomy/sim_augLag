#usr/bin/env python

### Testing smeared potentials in SIM suite.
### coded by KS

import time

import numpy as np

import sim


np.random.seed(12345)
NMol = 1

#define two atom types; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Charge = +1.0, Color = (1,0,0))
AtomTypeB = sim.chem.AtomType("B", Mass = 1.0, Charge = -1.0, Color = (0,1,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA, AtomTypeB])

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testsmearedcharge"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()


#set the system box length sizes
Sys.BoxL = [10.0, 10.0, 10.0]

#make a new potential energy term
P = sim.potential.SmearedCoulombEwCorr(Sys, Cut = 5.0, Shift = True, Coef=1, BornA=1, Filter=sim.atomselect.Pairs)
#P = sim.potential.SmearedCoulomb(Sys, Cut = 5.0, Shift = True, Coef=1, BornA=1, Filter=sim.atomselect.Pairs)

#add this term to the system forcefield
Sys.ForceField.append(P)
#NOTE: all potential energy terms must be added before compilation

#compile and load the system
Sys.Load()

#set initial positions and velocities
#sim.system.positions.CubicLattice(Sys, Random = 0.1)
rigidUnits = [x for x in Sys.RigidUnitsIter()]
N = len(rigidUnits)
Orient = np.zeros(Sys.Dim,float)
Sys.SetPos(rigidUnits[0],[0,0,0],Orient)
x = 2.5
Sys.SetPos(rigidUnits[1],[0,0,x],Orient)



Sys.Flags.CalcForce = True
Sys.Flags.CalcVirial = True
Sys.Flags.CalcDUParam = True
Sys.Flags.CalcDWParam = True
Sys.ForceField.Eval()
Sys.Flags.CalcsOff()
Sys.Flags.CalcVirial = True
PEnergy1 = Sys.PEnergy
Virial1 = Sys.Virial

P.SetTypeInd(1) #types are {0:AA, 1:AB, 2:BB}
print("Energy (Python): {}".format(Sys.ForceField[0].Val(x)))
print("Energy (Fortran):{}".format(Sys.PEnergy))

print("Virial (Python): {}".format(x * Sys.ForceField[0].DVal(x)))
print("Virial (Fortran):{}".format(Virial1))

'''
print "-", P.RPEnergy, P.SPEnergy, P.RVirial, Sys.PEnergy
Sys.TargetMol = 0
Sys.MolActive[0] = -1
Sys.ForceField.Eval(Mode = -2)
PEnergy2 = Sys.PEnergy
Virial2 = Sys.Virial
print "-", P.RPEnergy, P.SPEnergy, P.RVirial, Sys.PEnergy
Sys.ForceField.Eval()
PEnergy3 = Sys.PEnergy
Virial3 = Sys.Virial
print "-", P.RPEnergy, P.SPEnergy, P.RVirial, Sys.PEnergy
print PEnergy1, PEnergy2, PEnergy3
print Virial1, Virial2, Virial3
'''

