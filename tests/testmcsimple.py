#usr/bin/env python

### Testing for MC in SIM suite.
### coded by MSS

import numpy as np

import sim



#parameter for number of molecules
NMol = 250

#density
Rho = 0.80

#setpoint temperature
TempSet = 1.0

np.random.seed(12345)

#define two atom types; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.)

MolTypeA = sim.chem.MolType("MA", [AtomTypeA])
MolTypes = [MolTypeA]

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World(MolTypes, Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testmcsimple"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolTypeA.New()

#set the system box length sizes
Sys.BoxL = (len(Sys.Atom) / Rho)**(1./Sys.Dim)

print "Box volume is %.2f" % np.prod(Sys.BoxL)

#make a new potential energy term
Filter = sim.atomselect.InterPairs
P = sim.potential.LJ(Sys, Cut = 2.5, Filter = Filter, Epsilon = 1.0, Sigma = 1.0,  
                     Shift = True, Label = "LJ")
#add this term to the system forcefield
Sys.ForceField.append(P)

#NOTE: all potential energy terms must be added before compilation

#shortcut for the integrator  object
Int = Sys.Int

#compile and load the system
Sys.Load()

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys, Random = 0.1)
sim.system.velocities.Canonical(Sys, Temp = TempSet)
    
#run for some steps
Int.Method = Int.Methods.VVQuench
print "Minimizing"
Int.Run(2000)
Int.Reset()

print "Starting test" 

Sys.TempSet = TempSet

Int.Method = Int.Methods.MonteCarlo
Move = Int.Method.Moves[0]
Move.Delta = 0.1

#add verbose output of property averages
Sys.Measures.VerboseOutput(CycleFreq = 10)

print "CANONICAL:"
Int.Run(1000000) 
print "energy:", Sys.PEnergy
print "timing:", Int.TimeElapsed
print "average displacement acceptance:", Move.NAcc / float(Move.NAtt + 1.e-8)
print "average rotation acceptance:", Move.NAcc2 / float(Move.NAtt2 + 1.e-8)
print "\n"

Sys.MuSet[0] = -1
Move2 = Int.Method.Moves[1]
Move2.Weight = 1.

print "GRAND CANONICAL:"
Int.Run(1000000)
print "energy:", Sys.PEnergy
print "timing:", Int.TimeElapsed
print "average displacement acceptance:", Move.NAcc / float(Move.NAtt + 1.e-8)
print "average rotation acceptance:", Move.NAcc2 / float(Move.NAtt2 + 1.e-8)
print "average insert/delete acceptance:", Move2.NAcc / float(Move2.NAtt + 1.e-8)
print "\n"
Move2.Weight = 0.

Sys.PresSet = 1.
Move3 = Int.Method.Moves[2]
Move3.Weight = 1. / 50.
Move3.Delta = 10.

print "NPT:"
Int.Run(1000000)
print "energy:", Sys.PEnergy
print "timing:", Int.TimeElapsed
print "average displacement acceptance:", Move.NAcc / float(Move.NAtt + 1.e-8)
print "average rotation acceptance:", Move.NAcc2 / float(Move.NAtt2 + 1.e-8)
print "average volume change acceptance:", Move3.NAcc / float(Move3.NAtt + 1.e-8)
Move3.Weight = 0.



    