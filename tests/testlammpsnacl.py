#usr/bin/env python

### Testing for LAMMPS export using charges.
### coded by MSS

import numpy as np

import sim


#parameter for number of molecules
NMol = 100

#density  kg/m^3
Rho = 0.25 * 2160.
#convert to number density
Rho = Rho * 6.02e23 / 0.05844 / 1.e30

#setpoint temperature
TempSet = 1200.

np.random.seed(12345)

#define two atom types; give them names, masses, and charge
AtomTypeNa = sim.chem.AtomType("Na", Mass = 22.99, Charge = 1.) 
AtomTypeCl = sim.chem.AtomType("Cl", Mass = 35.45, Charge = -1.)

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("NaCl", [AtomTypeNa, AtomTypeCl])

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.AtomicUnits)

#make a system that exists in this world
SysName = "testlammpsnacl"
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
Cut = 5. * 2.5
P1 = sim.potential.LJCombine(Sys, Cut = Cut,
                             Filter = sim.atomselect.Pairs,
                             Epsilon = [2.77e-3, 0.1],
                             Sigma = [3.33045, 4.41724], 
                             Shift = True, Label = "LJ")
P2 = sim.potential.Ewald(Sys, Label = "EW", 
                         Cut = Cut, Shift = True) 
#add this term to the system forcefield
Sys.ForceField.extend([P1, P2])  
#NOTE: all potential energy terms must be added before compilation

#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate

#compile and load the system
Sys.Load()

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys)
sim.system.velocities.Canonical(Sys, Temp = TempSet)
Sys.TempSet = TempSet

#add verbose output of property averages
Sys.Measures.VerboseOutput(StepFreq = 1000)

#minimize energy with monte carlo
Int.Method = Int.Methods.MonteCarlo
print "Relaxing with Monte Carlo"
Int.Run(1000*NMol)
Int.Reset()

print "Running LAMMPS"

ret = sim.export.lammps.MakeLammpsTraj(Sys, Verbose = True,
                                       NStepsMin = 10000,
                                       NStepsEquil = 100000,
                                       NStepsProd = 100000,
                                       WriteFreq = 1000,
                                       Prefix = "testlmpnacl",
                                       DelTempFiles = False)
Traj, TrajFile = ret

print "frame", "SimPE", "LmpPE", "Err", "FracErr"
for (i, Pos) in enumerate(Traj):
    Sys.Pos = Pos
    Sys.ForceField.Eval()
    syspe = Sys.PEnergy
    trjpe = Traj.FrameData["PEnergy"]
    err = np.abs(syspe - trjpe)
    print i, syspe, trjpe, err, err/np.abs(syspe)