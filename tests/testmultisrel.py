#usr/bin/env python

### Testing for srel multi state minimization.
### coded by MSS

import numpy as np
import sim


TestNum = 1
UseLammps = False

StepsEquil = 10000
StepsProd = 100000
StepsStride = 100
np.random.seed(34231)


if UseLammps:
    OptClass = sim.srel.optimizetrajlammps.OptimizeTrajLammpsClass
else:
    OptClass = sim.srel.optimizetraj.OptimizeTrajClass

AtomType = sim.chem.AtomType("A", Mass = 1., Charge = 1.0)
MolType = sim.chem.MolType("M", [AtomType])
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

SysName = "multisrel%d" % TestNum
    
Sys = sim.system.System(World, Name = SysName)
for i in range(110):
    Sys += MolType.New()
    
BoxL1 = 125.**(1./3.)
BoxL2 = 138.**(1./3.)
Sys.BoxL[:] = BoxL1
Cut = 2.49999

if TestNum == 0:
    P = sim.potential.SoftSphere(Sys, Label = "SS", Cut = Cut,
                                 Filter = sim.atomselect.Pairs,
                                 Fixed = False, FixedExp = False,
                                 Exponent = 16, Sigma = 0.8, 
                                 Shift = True)
elif TestNum == 1:
    P = sim.potential.PairSpline(Sys, Label = "PS", Cut = Cut,
                                 Filter =sim.atomselect.Pairs,
                                 NKnot = 25)  
    #P.EmulateLJ(Sigma=1., Epsilon=1., Shift = True)
elif TestNum == 2:
    P = sim.potential.LJ(Sys, Label = "LJ", Cut = Cut, 
                         Filter = sim.atomselect.Pairs,
                         Epsilon = 1.1, Sigma = 1.05, 
                         Fixed = False, 
                         Shift = True)
else:
    raise ValueError("Not a recognized test number.")

Sys.ForceField.append(P)

#configure integrator
Int = Sys.Int
Int.Method.Thermostat = Int.Method.ThermostatLangevin


#set up the histograms
for P in Sys.ForceField:
    P.Arg.SetupHist(NBin = 10000, ReportNBin=100)

Sys.Load()

#initial positions and velocities
sim.system.positions.CubicLattice(Sys)
sim.system.velocities.Canonical(Sys, Temp = 1.0)
Sys.TempSet = 1.0  

#mapping function
Map = sim.atommap.PosMap()
for (i, a) in enumerate(Sys.Atom):
    Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]

#make two trajectories
Trj1 = sim.traj.Simple("../sampletraj/LJ-N110-V125-T10.trj.bz2")
Trj2 = sim.traj.Simple("../sampletraj/LJ-N110-V138-T10.trj.bz2")

Sys.ScaleBox(BoxL1)
Opt1 = OptClass(Sys, Map, Beta = 1., Traj = Trj1, TempFileDir = ".")
Sys.ScaleBox(BoxL2)
Opt2 = OptClass(Sys, Map, Beta = 1., Traj = Trj2, TempFileDir = ".")
OptList = [Opt1, Opt2]

for Opt in OptList:
    Opt.StepsEquil = StepsEquil
    Opt.StepsProd = StepsProd
    Opt.StepsStride = StepsStride

Optimizer = sim.srel.OptimizeMultiTrajClass(OptList)
Optimizer.FilePrefix = "%s" % Sys.Name
Optimizer.Run()





    