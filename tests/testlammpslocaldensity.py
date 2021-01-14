#usr/bin/env python

### Testing for local density potential in SIM suite.
### coded by MSS


import numpy as np

import sim


#parameter for number of molecules
NMol = 2

#box size
BoxL = 5.

#setpoint temperature
TempSet = 1.


np.random.seed(12345)


#define an atom type; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Charge = 0., Color = (1,0,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA])

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testmdlocaldensity"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = BoxL

#set the cutoff
Cut = 2.0
print "Cutoff is %f" % Cut

#make a new potential energy term
FilterAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA], Ordered = True)
P1 = sim.potential.LocalDensity(Sys, Cut = Cut, InnerCut = 0.5 * Cut,
                                Knots = [0., -1., -2., -3., -4., -5.], 
                                RhoMin = 0.0, RhoMax = 1.2, 
                                Label = "LD", Filter = FilterAA)
P2 = sim.potential.SoftSphere(Sys, Sigma = 1.0, Epsilon = 1.0, Label = "SS",
                              Cut = Cut, Shift = True, Filter = sim.atomselect.Pairs)
Sys.ForceField.extend([P1, P2])
#NOTE: all potential energy terms must be added before compilation

#set up the integrator 
Sys.TempSet = TempSet
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate
Int.Method.Thermostat = Int.Method.ThermostatLangevin

#define some things to measure.
#setup a measurement of a histogram of potential energy
Sys.Measures.PEnergy.SetupHist(0, 1000, 10)
#NOTE: all measurements must be defined (not necessarily active) before compilation

#compile and load the system
Sys.Load()

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys)
sim.system.velocities.Canonical(Sys, Temp = TempSet)


#print out local density
print "%11s %11s" % ("rij", "ene_locdens")
for rij in np.arange(0., Cut, 0.05):
    print "%11.5f %11.5f" % (rij, P1.Val(P1.RhoVal(rij)))


#test to export to lammps
Traj, TrajFile = sim.export.lammps.MakeLammpsTraj(Sys, Prefix = "testlmpld", 
                                                  NStepsEquil = 10000, NStepsProd = 10000, 
                                                  WriteFreq = 100)

 
print ("%3s " + "%11s "*8) % ("ind", "rij", "rho", "ene_locdens", "ene_total", "sim_ene", "lammps_ene", "diff", "sim_force")                                           
for (i, Pos) in enumerate(Traj):
    rij = sim.geom.Length(sim.geom.Minimage(Pos[0] - Pos[1], Sys.BoxL))
    SSEne = P2.Val(rij)
    rhoi = P1.RhoVal(rij)
    LDEne = P1.Val(rhoi)
    Sys.Pos = Pos
    Sys.Flags.CalcForce = True
    Sys.ForceField.Eval()
    SPE = Sys.PEnergy
    LPE = Traj.FrameData["PEnergy"]
    print ("%3d " + "%11.4e "*8) % (i, rij, rhoi, LDEne*2, SSEne+2*LDEne, SPE, LPE, abs(SPE-LPE), Sys.Atom[0].Force[0])
TrajEne, SysEne, TrajLogWeight, Nullx = sim.srel.base.GetEneTraj(Traj, Sys,
                                                                 ErrorNoEne = True, 
                                                                 ErrorDiffEne = True)


