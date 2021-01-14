#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS

import numpy as np
import sim
import time


sim.srel.base.DEBUG = False
sim.export.lammps.SkinDistance = 20.  

UseLammps = False

np.random.seed(34231)

#parameter for number of molecules
NMol = 1

#density
Rho = 0.001

#temp
TempSet = 1.0

#define one atom type
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Charge = 0.0, 
                              Color = (1,0,0), Radius = 0.3)
AtomTypeB = sim.chem.AtomType("B", Mass = 1.0, Charge = 0.0, 
                              Color = (0,1,0), Radius = 0.3)

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA, AtomTypeB] * 6)

#add a nonrigid bond between consecutive A atoms
for i in range(0, len(MolType) - 2, 2):
    MolType.Bond(i,i+2)
    
#add a nonrigid bond A and B atoms
for i in range(0, len(MolType), 2):
    MolType.Bond(i,i+1)

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testsrelpoly"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = (float(NMol*2) / Rho)**(1./Sys.Dim)
Sys.BoxL = 0.

#make new potential energy functions
FilterAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA], Bonded = True)
FilterAB = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeB], Bonded = True)
FilterAAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeA], Bonded = True)
FilterAAB = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeB], Bonded = True)
FilterAAAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeA, AtomTypeA], Bonded = True)
FilterAAAB = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeA, AtomTypeB], Bonded = True)
FilterNonBB = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB], MinBondOrd = 3)
PBondAA = sim.potential.Bond(Sys, Label = "BondAA", Filter = FilterAA, Dist0 = 1.0, FConst = 2)
PBondAB = sim.potential.Bond(Sys, Label = "BondAB", Filter = FilterAB, Dist0 = 0.5, FConst = 2)
PAngAAA  = sim.potential.AngleSpline(Sys, Label = "AngAAA", Filter = FilterAAA, NKnot = 20)
PAngAAB  = sim.potential.AngleSpline(Sys, Label = "AngAAB", Filter = FilterAAB, NKnot = 20)
PTorAAAA = sim.potential.TorsionSpline(Sys, Label = "TorAAAA", Filter = FilterAAAA, NKnot = 4)
PTorAAAB = sim.potential.TorsionSpline(Sys, Label = "TorAAAB", Filter = FilterAAAB, NKnot = 4)
PNonBB = sim.potential.PairSpline(Sys, Label = "NonBB", Filter = FilterNonBB, NKnot = 30, Cut = 2.5)
Sys.ForceField.extend([PBondAA, PBondAB, PAngAAA, PAngAAB, PTorAAAA, PTorAAAB, PNonBB])

#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate
Int.Method.Thermostat = Int.Method.ThermostatLangevin


#set up the histograms
for P in Sys.ForceField:
    P.Arg.SetupHist(NBin = 10000, ReportNBin=100)

Sys.Load()

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys, L = 10., Random = 0.1)
sim.system.velocities.Canonical(Sys, Temp = 1.0)

Map = sim.atommap.PosMap()
for (i, a) in enumerate(Sys.Atom):
    Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]

Trj = sim.traj.Simple("../sampletraj/polytraj.dat.bz2")

Sys.TempSet = 1.0

if UseLammps:
    OptClass = sim.srel.OptimizeTrajLammpsClass
    TempFileDir = "."  #for debug only
else:
    OptClass = sim.srel.OptimizeTrajClass
    TempFileDir = None
    
    
#sim.srel.base.DEBUG = True
Optimizer = OptClass(Sys, Map, Beta = 1., Traj = Trj, TempFileDir = TempFileDir)
Optimizer.FilePrefix = "%s_CG1" % Sys.Name
print Optimizer.Output0()

Optimizer.StepsEquil = 10000
Optimizer.StepsProd = 5000000   
Optimizer.StepsStride = 100


StartTime = time.time()

Stages = ["Bond", 
          "Bond AngleSpline PairSpline TorsionSpline"]


for Stage in Stages:
    PTypes = tuple([getattr(sim.potential, p) for p in Stage.split()])
    for P in Sys.ForceField:
        if isinstance(P, PTypes):
            P.UnfreezeParam()
            print "UNFREEZING", P.Name
        else:
            P.FreezeParam()
            print "FREEZING  ", P.Name
    Optimizer.Run()
  
        
print "\n\nTotal time: %.1f seconds" % (time.time() - StartTime)


    