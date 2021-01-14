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
TestNum = 1       #0 = soft sphere, 1=spline, 2=LJ
UseUPenalty = False
UseWPenalty = False
BlurSigma = 0.2
BlurNSample = 5


np.random.seed(34231)

#make system
AtomType = sim.chem.AtomType("A", Mass = 1., Charge = 1.0)
MolType = sim.chem.MolType("M", [AtomType])
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

SysName = "srel%d" % TestNum
    
Sys = sim.system.System(World, Name = SysName)
for i in range(110):
    Sys += MolType.New()
Sys.BoxL[:] = 5.

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
                         Epsilon = 1.0, Sigma = 1.0, 
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
Opt = sim.srel.OptimizeTrajClass(Sys, Map, Beta = 1., Traj = Trj, FilePrefix = "testsrel", 
                                 TempFileDir = TempFileDir, SaveLoadArgData = False,
                                 BlurSigma = BlurSigma, BlurNSample = BlurNSample)
if UseLammps:
    Opt = sim.srel.UseLammps(Opt)


#add penalties
if UseUPenalty:
    UPen = Opt.AddPenalty("PEnergy", -700., MeasureScale = 1./len(Sys.Atom))
if UseWPenalty:
    VPen = Opt.AddPenalty("Virial", -2400., MeasureScale = 1./len(Sys.Atom))


def ResetPotential():
    if TestNum == 0:
        P.SetParam(Sigma = 0.95, Exponent = 18.) 
    elif TestNum == 1:
        for (i,k) in enumerate(P.Knots):
            if not P.Knots.Fixed[i]:
                P.Knots[i] = 0
        if UseUPenalty or UseWPenalty:                
            Knots = [4.9894e+01 , 4.3644e+01 , 3.7394e+01 , 3.1144e+01 , 2.4894e+01 , 
                     1.8644e+01 , 1.2394e+01 , 6.1436e+00 , -1.1088e+00, -9.0677e-01, 
                     -7.8205e-01, -3.7126e-01, -2.7197e-01, -1.3825e-01, -8.8789e-02, 
                     -6.1139e-02, -3.7451e-02, -1.7963e-02, -1.2703e-02, -3.6489e-03]
            P.Knots = Knots    
            P.Update()
    elif TestNum == 2:
        P.SetParam(Sigma = 0.95, Epsilon = 4)

Opt.StepsEquil = 10000
Opt.StepsProd = 500000
Opt.StepsStride = 100
Opt.MinReweightFrames = None

StartTime = time.time()

   
if True:
    #conjugate gradient
    for Pen in Opt.Penalties:
        Pen.Coef = 0.
    Opt.Reset()
    ResetPotential()
    Opt.FilePrefix = Sys.Name
    Opt.Run()
    #Sys.ForceField.RelaxKnotConstraints()
    #Opt.Run()
                        

HasPenalty = UseUPenalty or UseWPenalty

if True and HasPenalty:
    #conjugate gradient
    Opt.Reset()
    ResetPotential()
    Opt.FilePrefix = Sys.Name
    Opt.RunStages(StageCoefs = [1.e-4, 1.e-2, 1.e-1, 1., 10., 100., 1000.]) 
    

print "\n\nTotal time: %.1f seconds" % (time.time() - StartTime)

