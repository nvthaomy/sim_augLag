#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS

import numpy as np

import sim
from sim.alg import bennett

#parameter for number of molecules
NMol = 110

#density
Rho = 0.8

#setpoint temperature
TempSet = 2.0

np.random.seed(12345)

#define two atom types; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.5, Charge = 0.0, Color = (1,0,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA])

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testbennett"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = (float(NMol) / Rho)**(1./Sys.Dim)

#make a new potential energy term
P = sim.potential.LJ(Sys, Cut = 2.5,
                     Filter = sim.atomselect.Pairs,
                     Epsilon = 1.0, Sigma = 1.0, 
                     Shift = True, Label = "LJ")
#add this term to the system forcefield
Sys.ForceField.append(P)

#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate

#compile and load the system
Sys.Load()

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys, Random = 0.1)
sim.system.velocities.Canonical(Sys, Temp = TempSet)

#run for some steps
Int.Method = Int.Methods.VVQuench
Int.Run(2000, ProgressText = "Minimizing")

Int.Method = Int.Methods.VVIntegrate
Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
Sys.TempSet = TempSet

Param1 = Sys.ForceField.Param.Val.copy()
P.SetParam(Epsilon = 1./1.9)
Param2 = Sys.ForceField.Param.Val.copy()
Sys.TempSet = 1.0

ret = bennett.CalcFEInt(Sys, 
                        StepsEquil = 50000, 
                        StepsProd = 500000, 
                        StepsStride = 10,
                        Param1 = Param1, Param2 = Param2,
                        MaxFEerr = 0.01,
                        CalcSrel = True)
NetFE, NetFEerr, Lambdas, DeltaFEs, DeltaFEerrs, Srel1, Srel2 = ret

print NetFE
print NetFEerr
print Lambdas
print DeltaFEs
print DeltaFEerrs

for (i, Lambda) in enumerate(Lambdas):
    ThisFE = sum(DeltaFEs[:i])
    print 1. * (1. - Lambda) + 1.9 * Lambda, ThisFE, Srel1[i], Srel2[i]