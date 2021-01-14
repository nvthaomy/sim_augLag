#usr/bin/env python

### Testing for LAMMPS in SIM suite.
### coded by MSS


import numpy as np

import sim

#parameter for number of molecules
NMol = 110

#density
Rho = 0.8

#setpoint temperature
TempSet = 1.0

np.random.seed(12345)

#define two atom types; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Color = (1,0,0))

#define a molecule type; give it a name and list of component atom types
MolTypeA = sim.chem.MolType("M", [AtomTypeA])

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolTypeA], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testlammps"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolTypeA.New()

#set the system box length sizes
Sys.BoxL = (float(NMol) / Rho)**(1./Sys.Dim)

#make a new potential energy term
P = sim.potential.LJ(Sys, Cut = 2.5,
                     Filter = sim.atomselect.Pairs,
                     Epsilon = 1.0, Sigma = 1.0, 
                     Shift = True, Label = "LJ")
#add this term to the system forcefield
Sys.ForceField.append(P)
#NOTE: all potential energy terms must be added before compilation

#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate

#define some things to measure.
#setup a measurement of a histogram of potential energy
Sys.Measures.PEnergy.SetupHist(-1200, -500, 300)
#NOTE: all measurements must be defined (not necessarily active) before compilation

#compile and load the system
Sys.Load()

#set initial positions and velocities
Sys.TempSet = TempSet
sim.system.positions.CubicLattice(Sys, Random = 0.1)
sim.system.velocities.Canonical(Sys, Temp = TempSet)

ret = sim.export.lammps.MakeLammpsTraj(Sys, Verbose = True,
                                       NStepsMin = 10000,
                                       NStepsEquil = 100000,
                                       NStepsProd = 500000,
                                       WriteFreq = 1000,
                                       Prefix = "testlmp")
Traj, TrajFile = ret

print "frame, sim PE, LAMMPS PE, absolute err, frac err"
for (i, Pos) in enumerate(Traj):
    Sys.Pos = Pos
    Sys.ForceField.Eval()
    syspe = Sys.PEnergy
    trjpe = Traj.FrameData["PEnergy"]
    err = np.abs(syspe - trjpe)
    print i, syspe, trjpe, err, err/np.abs(syspe)

