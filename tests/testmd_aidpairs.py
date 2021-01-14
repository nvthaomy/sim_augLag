#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS

import numpy as np
import sim

#parameter for number of molecules
NMol = 20

#density
Rho = 0.8

#setpoint temperature
TempSet = 2.0

#atomnames
AtomNames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']
AtomsPerMol = len(AtomNames)

np.random.seed(12345)

#define two atom types; give them names, masses, and charge
AtomTypes = []
for nm in AtomNames:
    at = sim.chem.AtomType(nm, Mass = 1.0)
    AtomTypes.append(at)

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", AtomTypes)

#add a rigid bond between component atoms #0 and 1, of length 1
#MolType.Bond(0,1,RLength = 1.)
#define the relative center of mass positions for the rigid molecule
#MolType.COMPos = np.array([[-0.5, 0, 0], [0.5, 0, 0]], float)
for i in range(AtomsPerMol-1):
    MolType.Bond(i, i+1)

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testmd"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = (float(NMol*AtomsPerMol) / Rho)**(1./Sys.Dim)

#make a new potential energy term
AIDPairs = np.zeros((AtomsPerMol, AtomsPerMol), dtype=bool)
for i in range(AtomsPerMol):
    AIDPairs[i,i] = True
    if i < AtomsPerMol-1: 
        AIDPairs[i,i+1] = True
        AIDPairs[i+1, i] = True
        
    
Filter = sim.atomselect.PolyFilter(Filters = [sim.atomselect.All, sim.atomselect.All], AIDPairs = AIDPairs)

Mode = 1
if Mode == 1:
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

sim.export.lammps.MakeLammpsMD(Sys, Prefix = "AIDPairs")


#add verbose output of property averages
Sys.Measures.VerboseOutput(StepFreq = 1000)

#run for some steps
Int.Method = Int.Methods.VVQuench
print "Minimizing"
Int.Run(2000)

NStep = 50000
Int.Method = Int.Methods.VVIntegrate
Sys.TempSet = TempSet

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


