#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS

import numpy as np
import sim
import time


sim.srel.base.DEBUG = False

#LAMMPS?
UseLammps = True


#TEST CONDITIONS
UseUPenalty = False
UseWPenalty = True


np.random.seed(34231)

NMol = 110

#make system
AtomType = sim.chem.AtomType("A", Mass = 1., Charge = 1.0)
MolType = sim.chem.MolType("M", [AtomType])
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

SysName = "srelld"
    
Sys = sim.system.System(World, Name = SysName)
for i in range(NMol):  
    Sys += MolType.New()
Sys.BoxL[:] = 5.

P = sim.potential.SoftSphere(Sys, Label = "SS", Cut = 2.5,
                             Filter = sim.atomselect.Pairs,
                             Fixed = True, FixedExp = True,
                             Exponent = 16, Sigma = 0.9, 
                             Shift = True)
LDFilter = sim.atomselect.PolyFilter([AtomType, AtomType], Ordered = True)
P2 = sim.potential.LocalDensity(Sys, Label = "LD", Cut = 1.5, InnerCut = 1.0,
                               Filter = LDFilter, NKnot = 25, RhoMin = 0., RhoMax = 50.)

Sys.ForceField.extend([P, P2])


#set up the histograms
for P in Sys.ForceField:
    P.Arg.SetupHist(NBin = 10000, ReportNBin=100)

Sys.Load()

#initial positions and velocities
sim.system.positions.CubicLattice(Sys)
sim.system.velocities.Canonical(Sys, Temp = 1.0)

#configure integrator
Int = Sys.Int

Int.Method = Int.Methods.VVQuench
Int.Run(1000, 'Minimizing')

Int.Method = Int.Methods.VVIntegrate        
Int.Method.Thermostat = Int.Method.ThermostatAndersen
Int.Method.AndersenStepFreq = 100
Int.Method.TimeStep = 0.001

Map = sim.atommap.PosMap()
for (i, a) in enumerate(Sys.Atom):
    Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]

Trj = sim.traj.Simple("../sampletraj/ljtraj.trj.bz2") 
Trj.BoxL = Sys.BoxL.copy()

Sys.TempSet = 1.0

#make srel Opt
TempFileDir = "."
Opt = sim.srel.OptimizeTrajClass(Sys, Map, Beta = 1., Traj = Trj, FilePrefix = "testsrelld", 
                                 TempFileDir = TempFileDir, SaveLoadArgData = False) 
if UseLammps:
    Opt = sim.srel.UseLammps(Opt)


#add penalties
if UseUPenalty:
    UPen = Opt.AddPenalty("PEnergy", -800., MeasureScale = 1./len(Sys.Atom))
if UseWPenalty:
    VPen = Opt.AddPenalty("Virial", -800., MeasureScale = 1./len(Sys.Atom))


Opt.StepsEquil = 10000
Opt.StepsProd = 250000
Opt.StepsStride = 100
Opt.MinReweightFrames = None

StartTime = time.time()

sim.srel.base.DEBUG = False 

#change to True to do a test of the objective function
if False:
    Opt.UpdateModTraj()  
    Opt.TestObj(DoDUParam = True, DoDWParam = True)
    exit()


#run first with no penalties
for Pen in Opt.Penalties:
    Pen.Coef = 0.
Opt.Reset()
Opt.FilePrefix = Sys.Name
Opt.Run()
    
                
HasPenalty = UseUPenalty or UseWPenalty
if True and HasPenalty:
    #conjugate gradient
    Opt.Reset()
    Opt.FilePrefix = Sys.Name
    Opt.RunStages(StageCoefs = [1.e-4, 1.e-2, 1.e-1, 1., 10., 100., 1000.]) 
    

print "\n\nTotal time: %.1f seconds" % (time.time() - StartTime)

