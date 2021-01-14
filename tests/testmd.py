#usr/bin/env python

### Testing for MD simulations in SIM suite.
### coded by MSS

import numpy as np

import sim

#parameter for number of molecules
NMol = 110

#density
Rho = 0.8799

#setpoint temperature and presure
TempSet = 2.0
PresSet = 8.0

#test
TestRigid = False
TestMode = 1

np.random.seed(12345)

#define two atom types; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.5, Charge = -0.6, Color = (1,0,0))
AtomTypeB = sim.chem.AtomType("B", Mass = 1.0, Charge = 0.6, Color = (0,1,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA, AtomTypeB])

if TestRigid:
    #add a rigid bond between component atoms #0 and 1, of length 1
    MolType.Bond(0,1,RLength = 1.)
    #define the relative center of mass positions for the rigid molecule
    MolType.COMPos = np.array([[-0.5, 0, 0], [0.5, 0, 0]], float)
else:    
    MolType.Bond(0,1)

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testmd"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = (float(NMol*2) / Rho)**(1./Sys.Dim)

#make a new potential energy term
Filter = sim.atomselect.Pairs
if TestMode == 1:
    P = sim.potential.LJ(Sys, Cut = 2.5,
                         Filter = Filter,
                         Epsilon = 1.0, Sigma = 1.0, 
                         Shift = True, Label = "LJ")
else:
    P = sim.potential.PairSpline(Sys, Cut = 2.5, NKnot = 50,
                                 Filter = Filter, Label = "Spline")
    P.EmulateLJ(1,1,Shift=True)
          
#add this term to the system forcefield
Sys.ForceField.append(P)
#NOTE: all potential energy terms must be added before compilation

if not TestRigid:
    PBond = sim.potential.Bond(Sys, Filter = sim.atomselect.BondPairs,
                               Dist0=0.8, FConst=3, Label = "Bond")
    Sys.ForceField.append(PBond)

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
Sys.Measures.Show()

#run for some steps
Int.Method = Int.Methods.VVQuench
print "Minimizing"
Int.Run(2000)

NStep = 10000
Int.Method = Int.Methods.VVIntegrate
Sys.TempSet = TempSet
Sys.PresSet = PresSet

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

Int.Method.Thermostat = Int.Method.ThermostatAndersenMassive
Int.Method.AndersenStepFreq = 1000
Sys.TempSet = TempSet
Sys.Measures.Reset()
print "MD-AndersenMassive"
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"

Int.Method.Thermostat = Int.Method.ThermostatAndersenCollision
Sys.TempSet = TempSet
Sys.Measures.Reset()
print "MD-AndersenCollision"
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"

Int.Method.Thermostat = Int.Method.ThermostatNoseHoover 
Sys.TempSet = TempSet
Sys.Measures.Reset()
print "MD-Nose Hoover"
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"

Int.Method.Thermostat = Int.Method.ThermostatNone
Sys.Measures.Reset()
print "MD-NVE"
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"

Int.Method.Thermostat = Int.Method.ThermostatNoseHoover 
Int.Method.Barostat = Int.Method.BarostatMonteCarlo
Sys.TempSet = TempSet
Sys.PresSet = PresSet
Sys.Measures.Reset()
print "MD-Nose Hoover + MC barostat"
Int.Run(NStep)
print "timing:", Int.TimeElapsed
print "\n"

