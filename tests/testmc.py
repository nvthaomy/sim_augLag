#usr/bin/env python

### Testing for Monte Carlo in SIM suite.
### coded by MSS

import numpy as np

import sim


#TestNum = 1  -- monatomic and LJ interactions
#TestNum = 2  -- monatomic and LJ+LD interactions
#TestNum = 3  -- rigid diatomic and LJ interactions
#TestNum = 4  -- rigid diatomic and LJ+LD interactions 
TestNum = 1


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

if TestNum == 1 or TestNum == 2:
    #define a molecule type; give it a name and list of component atom types
    MolTypeB = sim.chem.MolType("MB", [AtomTypeA])
    MolTypes.append(MolTypeB)
elif TestNum == 3 or TestNum == 4:
    #define a molecule type; give it a name and list of component atom types
    COMPos = np.array([[-0.25, 0, 0], [0.25, 0, 0]], dtype = float)
    MolTypeB = sim.chem.MolType("MB", [AtomTypeA, AtomTypeA], Rigid = True, COMPos = COMPos)
    MolTypes.append(MolTypeB)

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World(MolTypes, Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testmc"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
#for i in range(NMol):
#    Sys += MolType.New()
for i in range(NMol):
    if i % 2 == 0 or TestNum == 0:
        Sys += MolTypeA.New()
    else:
        Sys += MolTypeB.New()

#set the system box length sizes
Sys.BoxL = (len(Sys.Atom) / Rho)**(1./Sys.Dim)

print "Box volume is %.2f" % np.prod(Sys.BoxL)

#make a new potential energy term
Filter = sim.atomselect.InterPairs
P = sim.potential.LJ(Sys, Cut = 2.5, Filter = Filter, Epsilon = 1.0, Sigma = 1.0,  
                     Shift = True, Label = "LJ")
#add this term to the system forcefield
Sys.ForceField.append(P)
if TestNum == 2 or TestNum == 4:
    Cut = 2.5
    FilterAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA], Ordered = True)
    P2 = sim.potential.LocalDensity(Sys, Cut = Cut, InnerCut = 0.8 * Cut,
                                    Knots = [0., 1., 3., 7., 12., 18.], 
                                    RhoMin = 0.0, RhoMax = 50.0, 
                                    Label = "LD", Filter = FilterAA)
    Sys.ForceField.append(P2)
#NOTE: all potential energy terms must be added before compilation

#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate

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

Sys.MolActive[55:] = -1
Sys.UpdateActive()
Sys.MuSet[0] = -1
if TestNum > 0: Sys.MuSet[1] = -3
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

print "NPT:"
Int.Run(1000000)

print "energy:", Sys.PEnergy
print "timing:", Int.TimeElapsed
print "average displacement acceptance:", Move.NAcc / float(Move.NAtt + 1.e-8)
print "average rotation acceptance:", Move.NAcc2 / float(Move.NAtt2 + 1.e-8)
print "average volume change acceptance:", Move3.NAcc / float(Move3.NAtt + 1.e-8)
Move3.Weight = 0.

Sys.MuSet[0] = 0.
Sys.MuSet[1] = 0.
Move4 = Int.Method.Moves[3]
Move4.Weight = 1. / 50.
Move4.MID1 = 0
Move4.MID2 = 1
Sys.MolActive[:NMol/2] = 1
Sys.MolActive[NMol/2:] = -1
Sys.UpdateActive()

print "SEMI GRAND:"
Int.Run(1000000)
print "energy:", Sys.PEnergy
print "timing:", Int.TimeElapsed
print "average displacement acceptance:", Move.NAcc / float(Move.NAtt + 1.e-8)
print "average rotation acceptance:", Move.NAcc2 / float(Move.NAtt2 + 1.e-8)
print "average mutation acceptance:", Move4.NAcc / float(Move4.NAtt + 1.e-8)

print "SEMI GRAND WITH 10 MULTICANONICAL WEIGHT UPDATES:"
Move4.ResetTM()
for i in range(10):
    Int.Run(1000000)
    ThisWeights = Move4.CalcTMWeights()
    Move4.BoltzWeights = -ThisWeights
 
print "State i, Hist[i], Weights[i]"
for (i, iTM) in enumerate(Move4.TM):
    print i, np.sum(iTM), ThisWeights[i]
    