#usr/bin/env python

### Testing for SPC/E model in SIM suite.
### coded by MSS

import numpy as np

import sim

#parameter for number of molecules
NMol = 150

#density  kg/m^3
Rho = 1000.
#convert to number density
Rho = Rho * 6.02e23 / 0.018 / 1.e30

#setpoint temperature
TempSet = 298.

np.random.seed(12345)

#define two atom types; give them names, masses, and charge
AtomTypeO = sim.chem.AtomType("O", Mass = 15.99, Charge = -0.8476, Color = (1,0,0))
AtomTypeH = sim.chem.AtomType("H", Mass = 1.00794, Charge = 0.4238, Color = (0,1,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("W", [AtomTypeO, AtomTypeH, AtomTypeH], Rigid = True)

#define the relative center of mass positions for the rigid molecule
Pos = np.array([[0.00000,  0.00000,  0.00000],
                [0.81649,  0.57736,  0.00000], 
                [-0.81649,  0.57736,  0.00000]], dtype=float)
#add rigid bonds
bl01 = sim.geom.Length(Pos[0] - Pos[1])
bl02 = sim.geom.Length(Pos[0] - Pos[2])
bl12 = sim.geom.Length(Pos[1] - Pos[2])
MolType.Bond(0,1,RLength = bl01)
MolType.Bond(0,2,RLength = bl02)
MolType.Bond(1,2,RLength = bl12)
MolType.COMPos = Pos


#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.AtomicUnits)

#make a system that exists in this world
SysName = "testmdspceatomic"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
L = (float(NMol) / Rho)**(1./3)
Sys.BoxL = L
print "Box dimensions: ", Sys.BoxL
print "Box volume: ", np.prod(Sys.BoxL)

#make a new potential energy term
Filter = sim.atomselect.PolyFilter([AtomTypeO, AtomTypeO])
P1 = sim.potential.LJ(Sys, Cut = 3.1658 * 2.5,
                      Filter = Filter,
                      Epsilon = 0.65, Sigma = 3.1658, 
                      Shift = True, Label = "LJ")
P2 = sim.potential.Ewald(Sys, Label = "EW", ExcludeBondOrd = 3, 
                         Cut = 3.1658 * 2.5, Shift = True)
#add this term to the system forcefield
Sys.ForceField.extend([P1, P2])
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
sim.system.positions.CubicLattice(Sys)
sim.system.velocities.Canonical(Sys, Temp = TempSet)

#add verbose output of property averages
Sys.Measures.VerboseOutput(StepFreq = 1000)

#minimize energy with monte carlo
Int.Method = Int.Methods.MonteCarlo
print "Relaxing with Monte Carlo"
Int.Run(1000*NMol)
Int.Reset()

NStep = 10000
Int.Method = Int.Methods.VVIntegrate
#set timestep
Int.Method.TimeStep = 0.01
Int.Method.Thermostat = Int.Method.ThermostatLangevin
Sys.TempSet = TempSet
Sys.Measures.Reset()
Int.Reset()
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

Int.Method = Int.Methods.MonteCarlo
Int.Reset()
Sys.Measures.Reset()
print "MC"
Int.Run(NStep*NMol)
print "timing:", Int.TimeElapsed
print "\n"

