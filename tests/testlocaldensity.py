#usr/bin/env python

### Testing for local density potential in SIM suite.
### coded by MSS

import time

import numpy as np

import sim

#parameter for number of molecules
NMol = 100

#density
Rho = 0.4

#setpoint temperature
TempSet = 1.

np.random.seed(12345)

#define an atom type; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Charge = 0., Color = (1,0,0))
AtomTypeB = sim.chem.AtomType("B", Mass = 1.0, Charge = 0., Color = (0,1,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA,AtomTypeB])

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testmdlocaldensity"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = (float(NMol) / Rho)**(1./Sys.Dim)

#set the cutoff
Cut = Sys.BoxL[0] * 0.5
print "Cutoff is %f" % Cut

#make a new potential energy term
FilterAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA], Ordered = True)
P1 = sim.potential.LocalDensity(Sys, Cut = Cut, InnerCut = 0.5 * Cut,
                                Knots = [0., 1., 3., 7., 12., 18.], 
                                RhoMin = 0.0, RhoMax = 50., 
                                Label = "LD", Filter = FilterAA)
P2 = sim.potential.LJ(Sys, Sigma = 1.0, Epsilon = 1.0, Label = "LJ",
                      Cut = Cut, Shift = True, Filter = sim.atomselect.Pairs)
Sys.ForceField.extend([P1, P2])
#NOTE: all potential energy terms must be added before compilation

#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate

#define some things to measure.
#setup a measurement of a histogram of potential energy
Sys.Measures.PEnergy.SetupHist(0, 1000, 10)
#NOTE: all measurements must be defined (not necessarily active) before compilation

#compile and load the system
Sys.Load()

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys)
sim.system.velocities.Canonical(Sys, Temp = TempSet)

#add an action to happen during integration
def PrintMeasures(Sys):
    print Sys.Measures.ReportVals()
Act = sim.integrate.Action(StepFreq = 100, Fn = PrintMeasures)
Int.Actions.append(Act)


#run for some steps
Int.Method = Int.Methods.VVQuench
print "Minimizing"
print Sys.Measures.ReportHead()
Int.Run(1000)

NStep = 5000

Int.Method = Int.Methods.VVIntegrate
Int.Method.Thermostat = Int.Method.ThermostatAndersen
Int.Method.AndersenStepFreq = 200
Sys.TempSet = TempSet
Sys.Measures.Reset()
t1 = time.time()
print "MD-Andersen"
print Sys.Measures.ReportHead()
Int.Run(NStep)
print "-"*len(Sys.Measures.ReportHead())
Sys.Measures.Reset()
Int.Run(NStep)
print Sys.Measures.ReportAvgs()
Sys.Measures.Reset()

print "-"*len(Sys.Measures.ReportHead())
Int.Method.Thermostat = Int.Method.ThermostatNone
print "MD-NVE"
print Sys.Measures.ReportHead()
Int.Run(NStep)
print Sys.Measures.ReportAvgs()
#print "\n", Sys.Measures.PEnergy.HistString()
t2 = time.time()

print "timing:", t2-t1


