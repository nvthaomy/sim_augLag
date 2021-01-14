#usr/bin/env python

### Testing for smoothed torsion potential in SIM suite.
### coded by MSS

import time

import numpy as np

import sim


#parameter for number of molecules
NMol = 1

#density
Rho = 0.001

#temp
TempSet = 1.0

np.random.seed(12345)

#define one atom type
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Charge = 0.0, 
                              Color = (1,0,0), Radius = 0.3)

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA] * 4)

#add a nonrigid bond between consecutive A atoms
for i in range(0, len(MolType) - 1):
    MolType.Bond(i,i+1)
    
#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3)

#make a system that exists in this world
SysName = "testmdsmoothedtorsion"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = (float(NMol) / Rho)**(1./Sys.Dim)
#Sys.BoxL = 0.

#make new potential energy functions
FilterAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA], Bonded = True)
FilterAAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeA], Bonded = True)
FilterAAAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeA, AtomTypeA], Bonded = True)
PBondAA = sim.potential.Bond(Sys, Filter = FilterAA, Dist0 = 1.0, FConst = 2)
PTorAAAA = sim.potential.SmoothedTorsionSpline(Sys, Filter = FilterAAAA,
                                       Knots = [1., .5, 0., 0.5])
Sys.ForceField.extend([PBondAA, PTorAAAA])
                        
#NOTE: all potential energy terms must be added before compilation

#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate
Int.Method.TimeStep = 0.001

#set up the histograms
for P in Sys.ForceField:
    P.Arg.SetupHist(NBin = 10000, ReportNBin=100) 

#compile and load the system
Sys.Lib.ForceCompile = True
Sys.Load()

#set initial positions and velocities
sim.system.init.positions.CubicLattice(Sys, L = 10., Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = 1.0)

#report measurement values
Sys.Measures.VerboseOutput(StepFreq = 10000)

#set visualizer and add an update action to the integrator
#import sim.system.visualize
#Vis = sim.system.visualize.Visualizer3D(Sys, ShowBonds = True, Label = True)
#Vis.AddAction(Int, TimeFreq = 0.2)

#run for some steps
Int.Method = Int.Methods.VVQuench
print "Minimizing"
Int.Run(10000)

NStep = 100000
NStepEquil = 100000
Int.Method = Int.Methods.VVIntegrate
Sys.TempSet = TempSet

#FOR PROFILING PURPOSES
#t1 = time.time()
#Int.Method.Thermostat = Int.Method.ThermostatLangevin
#Int.Method.LangevinGamma = 10
#Sys.TempSet = TempSet
#Sys.Measures.Reset()
#print "MD-Langevin-Profile"
#Int.Run(50000)
#print time.time() - t1
#exit()

Int.Method.Thermostat = Int.Method.ThermostatLangevin
Int.Method.LangevinGamma = 10
Sys.TempSet = TempSet
Sys.Measures.Reset()
print "MD-Langevin"
Int.Run(NStepEquil)
Sys.Measures.Reset()
t1 = time.time()
Int.Run(NStep)
print "timing:", time.time()-t1
print "\n"
Sys.Measures.Reset()


#run NVE and make trajectory to save
trj = sim.traj.raw.RawWrite("testmdpoly.traj.dat")
trj.AddAction(Int, StepFreq = 50000)
Int.Method.Thermostat = Int.Method.ThermostatNone
Sys.Measures.Reset()
t1 = time.time()
print "MD-NVE"
Int.Run(NStep)
print "timing:", time.time()-t1
print "\n"
trj.DelAction(Int)
trj.Close()

trj = sim.traj.raw.Raw("testmdpoly.traj.dat")
print "Length of trajectory: %d\n" % len(trj)
Beta = 1.0

#reparse trajectory to get averages
print "Reparsing trajectory for averages...\n"
ret = sim.srel.base.CalcAvgsTraj(trj, Sys, Beta, Reweight = False, CalcDeriv = True)
AvgDUParam1, AvgSqDUParam, AvgDDUParam1, Count, ReweightFrac, Srel2Term, AvgPEnergy1 = ret

trj.Reset()

#reparse trajectory to get histograms
print "Reparsing trajectory for histograms...\n"
RefPos, RefBoxL = sim.srel.base.ParseArgTraj(trj, Sys, Beta)

#recalculate averages
print "Reparsing histograms for averages...\n"
AvgDUParam2, AvgDDUParam2, AvgPEnergy2 = sim.srel.base.CalcAvgsHist(Sys)

#print "Bond distance histogram:"
#print PBondAA.HistString()
#print "\n"

print "Comparison of penergy averages:"
print "%11.4e %11.4e" % (AvgPEnergy1, AvgPEnergy2)

print "Comparison of duparam averages:"
for (i, (x,y)) in enumerate(zip(AvgDUParam1,AvgDUParam2)):
    print "%3d %11.4e %11.4e %11.4e" % (i, x, y, abs(x-y))
    