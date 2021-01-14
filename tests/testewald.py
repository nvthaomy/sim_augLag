#usr/bin/env python

### Testing for Ewald in SIM suite.
### coded by MSS

import numpy as np

import sim


np.random.seed(12345)
NMol = 2

#define two atom types; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Charge = +0.72, Color = (1,0,0))
AtomTypeB = sim.chem.AtomType("B", Mass = 1.0, Charge = -0.72, Color = (0,1,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA, AtomTypeB])

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.AtomicUnits)

#make a system that exists in this world
SysName = "testewald"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = [2.15443469003188, 2.15443469003188, 2.15443469003188]

#make a new potential energy term
P = sim.potential.Ewald(Sys, Cut = 2.5, Shift = True)
#add this term to the system forcefield
Sys.ForceField.append(P)
#NOTE: all potential energy terms must be added before compilation

#compile and load the system
Sys.Load()

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys, Random = 0.1)

Sys.Flags.CalcForce = True
Sys.Flags.CalcVirial = True
Sys.Flags.CalcDUParam = True
Sys.Flags.CalcDWParam = True
Sys.ForceField.Eval()
Sys.Flags.CalcsOff()
PEnergy1 = Sys.PEnergy
print "-", P.RPEnergy, P.SPEnergy, Sys.PEnergy
Sys.TargetMol = 0
Sys.MolActive[0] = -1
Sys.ForceField.Eval(Mode = -2)
PEnergy2 = Sys.PEnergy
print "-", P.RPEnergy, P.SPEnergy, Sys.PEnergy
Sys.ForceField.Eval()
PEnergy3 = Sys.PEnergy
print "-", P.RPEnergy, P.SPEnergy, Sys.PEnergy
print PEnergy1, PEnergy2, PEnergy3


