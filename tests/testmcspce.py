#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS

import time

import numpy as np

import sim


#parameter for number of molecules
NMol = 150

#density  kg/m^3
Rho = 1000.
#convert to number density
Rho = Rho * 6.02e23 / 0.018 / 1.e30

#setpoint temperature
TempSet = 298.

np.random.seed(12345)

#define two atom types; give them names, masses, and charge
AtomTypeO = sim.chem.AtomType("O", Mass = 15.99, Charge = -0.8476, Color = (1,0,0))
AtomTypeH = sim.chem.AtomType("H", Mass = 1.00794, Charge = 0.4238, Color = (0,1,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("W", [AtomTypeO, AtomTypeH, AtomTypeH], Rigid = True)

#define the relative center of mass positions for the rigid molecule
Pos = np.array([[0.00000,  0.00000,  0.00000],
                [0.81649,  0.57736,  0.00000], 
                [-0.81649,  0.57736,  0.00000]], dtype=float)
#add rigid bonds
bl01 = sim.geom.Length(Pos[0] - Pos[1])
bl02 = sim.geom.Length(Pos[0] - Pos[2])
bl12 = sim.geom.Length(Pos[1] - Pos[2])
MolType.Bond(0,1,RLength = bl01)
MolType.Bond(0,2,RLength = bl02)
MolType.Bond(1,2,RLength = bl12)
MolType.COMPos = Pos


#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.AtomicUnits)

#make a system that exists in this world
SysName = "testmcspce"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
L = (float(NMol) / Rho)**(1./3)
Sys.BoxL = L
print "Box dimensions: ", Sys.BoxL
print "Box volume: ", np.prod(Sys.BoxL)

#make a new potential energy term
Filter = sim.atomselect.PolyFilter([AtomTypeO, AtomTypeO])
P1 = sim.potential.LJ(Sys, Cut = 3.1658 * 2.5,
                      Filter = Filter,
                      Epsilon = 0.65, Sigma = 3.1658, 
                      Shift = True, Label = "LJ")
P2 = sim.potential.Ewald(Sys, Label = "EW", ExcludeBondOrd = 3, 
                         Cut = 3.1658 * 2.5, Shift = True)
#add this term to the system forcefield
Sys.ForceField.extend([P1, P2])
#NOTE: all potential energy terms must be added before compilation


#compile and load the system
Sys.Load()

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys)
sim.system.velocities.Canonical(Sys, Temp = TempSet)
    
print "Starting test" 

Sys.TempSet = TempSet


print "ENERGY MINIMIZING (BRIEF)"
#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVQuench
Int.Run(1000)


NStep = 1000000

Int.Method = Int.Methods.MonteCarlo
Move = Int.Method.Moves[0]
Move.Delta = 1.
Move.Delta = 0.1
Move.Delta2 = np.pi / 18.

print "FINDING MAX DISPLACEMENTS"
for i in range(10):
    Int.Run(10000)
    print "Current acceptance fractions: ", Move.AccFrac()
    print "   current max displacements: ", Move.Delta, Move.Delta2
    Move.UpdateDelta()
    print "   updated max displacements: ", Move.Delta, Move.Delta2
    print ""
    Int.Reset()

#add verbose output of property averages
Sys.Measures.VerboseOutput(CycleFreq = 1)

print "CANONICAL:"
Int.Run(NStep) 
print "energy:", Sys.PEnergy
print "timing:", Int.TimeElapsed
acc1, acc2 = Move.AccFrac()
print "average displacement acceptance:", acc1
print "average rotation acceptance:", acc2
print "\n"


Sys.MuSet[0] = -1
Move2 = Int.Method.Moves[1]
Move2.Weight = 1.

print "GRAND CANONICAL:"
Int.Run(NStep)
print "energy:", Sys.PEnergy
print "timing:", Int.TimeElapsed
acc1, acc2 = Move.AccFrac()
acc3 = Move2.AccFrac()[0]
print "average displacement acceptance:", acc1
print "average rotation acceptance:", acc2
print "average insert/delete acceptance:", acc3
print "\n"
Move2.Weight = 0.


Sys.PresSet = 0.
Move3 = Int.Method.Moves[2]
Move3.Delta = 3.
Move3.Weight = 1. / 50.

print "NPT:"
Int.Run(NStep)
print "energy:", Sys.PEnergy
print "timing:", Int.TimeElapsed
acc1, acc2 = Move.AccFrac()
acc3 = Move3.AccFrac()[0]
print "average displacement acceptance:", acc1
print "average rotation acceptance:", acc2
print "average volume change acceptance:", acc3
Move3.Weight = 0.
