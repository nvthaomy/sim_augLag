#usr/bin/env python

### Testing for MD for LJ system in SIM suite.
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
SysName = "testmdlj"
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
sim.system.positions.CubicLattice(Sys, Random = 0.1)
sim.system.velocities.Canonical(Sys, Temp = TempSet)

#add verbose output of property averages
Sys.Measures.VerboseOutput(StepFreq = 1000)

#run for some steps
Int.Method = Int.Methods.VVQuench
print "Minimizing"
Int.Run(2000)

NStep = 10000
Int.Method = Int.Methods.VVIntegrate
Sys.TempSet = TempSet

print "Finding an appropriate time step...."
Sys.Measures.VerboseOutput(StepFreq = 0)
sim.integrate.velverlet.FindTimeStep(Sys, NSteps = 10000)
Sys.Measures.VerboseOutput(StepFreq = 1000)
Sys.Measures.Reset()
Sys.Int.Reset()

print "Now running production runs...."
Int.Method.Thermostat = Int.Method.ThermostatLangevin
Int.Method.LangevinGamma = 10
Sys.TempSet = TempSet
Sys.Measures.Reset()
print "MD-Langevin"
Int.Run(NStep)
Sys.Measures.Reset()
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"
Sys.Measures.Reset()

Int.Method.Thermostat = Int.Method.ThermostatAndersenMassive
Int.Method.AndersenStepFreq = 1000
Sys.TempSet = TempSet
Sys.Measures.Reset()
print "MD-AndersenMassive"
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"
Sys.Measures.Reset()

Int.Method.Thermostat = Int.Method.ThermostatAndersenCollision
Sys.TempSet = TempSet
Sys.Measures.Reset()
print "MD-AndersenCollision"
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"
Sys.Measures.Reset()

Int.Method.Thermostat = Int.Method.ThermostatNoseHoover 
Sys.TempSet = TempSet
Sys.Measures.Reset()
print "MD-Nose Hoover"
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"
Sys.Measures.Reset()

Int.Method.Thermostat = Int.Method.ThermostatNone
Sys.Measures.Reset()
print "MD-NVE"
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"


