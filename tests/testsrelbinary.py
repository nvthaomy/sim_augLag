#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS


import numpy as np
import sim


np.random.seed(3423)

AtomTypeA = sim.chem.AtomType("A")
AtomTypeB = sim.chem.AtomType("B")
MolTypeA = sim.chem.MolType("A", [AtomTypeA])
MolTypeB = sim.chem.MolType("B", [AtomTypeB])
World = sim.chem.World([MolTypeA, MolTypeB], Dim = 3, Units = sim.units.DimensionlessUnits)

SysName = "srelbinary"
    
Sys = sim.system.System(World, Name = SysName)
for i in range(250):
    Sys += MolTypeA.New()
for i in range(250):
    Sys += MolTypeB.New()
Sys.BoxL[:] = 6.9336

PAA = sim.potential.SoftSphere(Sys, Label = "AA", Cut = 2.5,
                               Filter = sim.atomselect.PolyFilter([MolTypeA, MolTypeA]),
                               Fixed = False, FixedExp = False,
                               Exponent = 2, Sigma = 1, 
                               Shift = True)
PAB = sim.potential.SoftSphere(Sys, Label = "AB", Cut = 2.5,
                               Filter =  sim.atomselect.PolyFilter([MolTypeA, MolTypeB]),
                               Fixed = False, FixedExp = False,
                               Exponent = 2, Sigma = 1, 
                               Shift = True)
PBB = sim.potential.SoftSphere(Sys, Label = "BB", Cut = 2.5,
                               Filter = sim.atomselect.PolyFilter([MolTypeB, MolTypeB]),
                               Fixed = False, FixedExp = False,
                               Exponent = 2, Sigma = 1, 
                               Shift = True)


Sys.ForceField.append(PAA)
Sys.ForceField.append(PAB)
Sys.ForceField.append(PBB)
#set up the histograms
for P in Sys.ForceField:
    P.Arg.SetupHist(NBin = 10000, ReportNBin=100)
    
Sys.Load()

Map = sim.atommap.PosMap()
for (i, a) in enumerate(Sys.Atom):
    Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]

Trj = sim.traj.Simple("../sampletraj/m12-r1.5-T1.00.trj.bz2")

Opt = sim.srel.OptimizeTrajClass(Sys, Map, Beta = 1., Traj = Trj, FilePrefix = "srelbinary")

#setup the run
Opt.StepsEquil = 5000
Opt.StepsProd = 10000 
Opt.StepsStride = 100
Opt.MinReweightFrames = 50
Int = Sys.Int
Int.Method = Int.Methods.VVIntegrate        
Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
Opt.Run()




    