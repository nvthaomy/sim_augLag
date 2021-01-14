#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS

import numpy as np
import sim
import time


sim.srel.base.DEBUG = False

#integrator
IntMode = 0  #0 for MD, 1 for MC

#TEST CONDITIONS
TestNum = 1        #0 = soft sphere, 1=spline, 2=LJ

#parameter for number of molecules
NMol = 128

#density
Rho = 0.80

#setpoint temperature and pressure
TempSet = 1.0
PresSet = 2.0


np.random.seed(34231)

#make system
AtomType = sim.chem.AtomType("A", Mass = 1., Charge = 1.0)
MolType = sim.chem.MolType("M", [AtomType])
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

SysName = "srelnpt%d" % TestNum
    
Sys = sim.system.System(World, Name = SysName)
for i in range(NMol):
    Sys += MolType.New()
#set the system box length sizes
Sys.BoxL = (len(Sys.Atom) / Rho)**(1./Sys.Dim)

if TestNum == 0:
    P = sim.potential.SoftSphere(Sys, Label = "SS", Cut = 2.5,
                                 Filter = sim.atomselect.Pairs,
                                 Fixed = False, FixedExp = False,
                                 Exponent = 16, Sigma = 0.9, 
                                 Shift = True)
elif TestNum == 1:
    P = sim.potential.PairSpline(Sys, Label = "PS", Cut = 2.5,
                                 Filter = sim.atomselect.Pairs,
                                 NKnot = 20)  
elif TestNum == 2:
    P = sim.potential.LJ(Sys, Label = "LJ", Cut = 2.5, 
                         Filter = sim.atomselect.Pairs,
                         Epsilon = 0.8, Sigma = 0.9, 
                         Shift = True)
else:
    raise ValueError("Not a recognized test number.")

Sys.ForceField.append(P)


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
Int.Reset()

Sys.TempSet = TempSet 
Sys.PresSet = PresSet  

if IntMode == 0:
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
    Int.Method.Barostat = Int.Method.BarostatMonteCarlo
    StepsScale = 1
else:
    Int.Method = Int.Methods.MonteCarlo   
    MoveDisplace = Int.Method.Moves[0]
    MoveDisplace.Delta = 0.1
    MoveV = Int.Method.Moves[2]
    MoveV.Weight = 1./50.
    StepsScale = NMol

Map = sim.atommap.PosMap()
for (i, a) in enumerate(Sys.Atom):
    Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]

Trj = sim.traj.Lammps("../sampletraj/LJ-N128-P2-T1.lammpstrj.bz2")  

#make srel Opt
TempFileDir = "."
Opt = sim.srel.OptimizeTrajClass(Sys, Map, Beta = 1., Traj = Trj, 
                                 FilePrefix = "testsrelnpt%d" % TestNum, 
                                 TempFileDir = TempFileDir)  

Opt.StepsEquil = 10000 * StepsScale
Opt.StepsProd = 100000 * StepsScale
Opt.StepsStride = 100 * StepsScale

StartTime = time.time()

Opt.Run()    

print "\n\nTotal time: %.1f seconds" % (time.time() - StartTime)

