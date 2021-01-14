#usr/bin/env python

### Testing for a model polymer with full bonded interactions in SIM suite.
### coded by MSS

import numpy as np

import sim

#write trajectory?
WriteTraj = False

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
AtomTypeB = sim.chem.AtomType("B", Mass = 1.0, Charge = 0.0, 
                              Color = (0,1,0), Radius = 0.3)

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA, AtomTypeB] * 6)

#add a nonrigid bond between consecutive A atoms
for i in range(0, len(MolType) - 2, 2):
    MolType.Bond(i,i+2)
    
#add a nonrigid bond A and B atoms
for i in range(0, len(MolType), 2):
    MolType.Bond(i,i+1)

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testmdpoly"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = (float(NMol*2) / Rho)**(1./Sys.Dim)

#make new potential energy functions
FilterA = sim.atomselect.PolyFilter([AtomTypeA])
FilterAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA], Bonded = True)
FilterAB = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeB], Bonded = True)
FilterAAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeA], Bonded = True)
FilterAAB = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeB], Bonded = True)
FilterAAAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeA, AtomTypeA], Bonded = True)
FilterAAAB = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA, AtomTypeA, AtomTypeB], Bonded = True)
FilterNonBB = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB], MinBondOrd = 3)
PNonBB = sim.potential.LJ(Sys, Filter = FilterNonBB, Epsilon = 1., Sigma = 1., Cut = 4.)
PBondAA = sim.potential.Bond(Sys, Filter = FilterAA, Dist0 = 1.0, FConst = 2)
PBondAB = sim.potential.Bond(Sys, Filter = FilterAB, Dist0 = 0.5, FConst = 2)
Sys.ForceField.extend([PBondAA, PBondAB, PNonBB])
PAngAAA  = sim.potential.Angle(Sys, Filter = FilterAAA,
                               Theta0 = 120. * np.pi / 180., FConst = 2.0)
PAngAAB  = sim.potential.Angle(Sys, Filter = FilterAAB,
                               Theta0 = 90. * np.pi / 180., FConst = 2.0)
PTorAAAA = sim.potential.TorsionSpline(Sys, Filter = FilterAAAA,
                                       Knots = [1., .5, 0., 0.5])
PTorAAAB = sim.potential.TorsionSpline(Sys, Filter = FilterAAAB,
                                       Knots = [0., 1.0, 0., 0.5])
Sys.ForceField.extend([PAngAAA, PAngAAB, PTorAAAA, PTorAAAB])

                        
#NOTE: all potential energy terms must be added before compilation

#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate
Int.Method.TimeStep = 0.001

Filters = [sim.atomselect.Filter(Types = MolType[i]) for i in range(0,8,2)]
Filter = sim.atomselect.PolyFilter(Filters = Filters, Bonded = True)
Sys.Measures.append(sim.measure.Dihedral(Sys, Filter, StepFreq = 10))
Sys.Measures[-1].Report = False
Filters = [sim.atomselect.Filter(Types = MolType[i]) for i in range(0,6,2)]
Filter = sim.atomselect.PolyFilter(Filters = Filters, Bonded = True)
Sys.Measures.append(sim.measure.Angle(Sys, Filter, StepFreq = 10))
Sys.Measures[-1].Report = False
Filters = [sim.atomselect.Filter(Types = MolType[i]) for i in range(0,4,2)]
Filter = sim.atomselect.PolyFilter(Filters = Filters, Bonded = True)
Sys.Measures.append(sim.measure.Distance(Sys, Filter, StepFreq = 10))
Sys.Measures[-1].Report = False

#test filter
for P in Sys.ForceField:
    print "Potential %s types:" % str(P)
    for STypes in P.Filter.Select(Sys.World.SiteTypes):
        print "  ", STypes
print "\n"


#compile and load the system
Sys.Lib.ForceCompile = True
Sys.Load()

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys, L = 10., Random = 0.1)
sim.system.velocities.Canonical(Sys, Temp = 1.0)

#report measurement values
Sys.Measures.VerboseOutput(StepFreq = 10000)


#run for some steps
Int.Method = Int.Methods.VVQuench
print "Minimizing"
Int.Run(10000)

NStep = 100000
NStepEquil = 100000
Int.Method = Int.Methods.VVIntegrate
Sys.TempSet = TempSet

Int.Method.Thermostat = Int.Method.ThermostatLangevin
Int.Method.LangevinGamma = 10
Sys.TempSet = TempSet
Sys.Measures.Reset()
print "MD-Langevin"
Int.Run(NStepEquil)
Sys.Measures.Reset()

if WriteTraj:
    Trj = sim.traj.simple.SimpleWrite("polytraj.dat.bz2")
    Trj.AddAction(Int, StepFreq = 100)
Int.Run(NStep)
if WriteTraj:
    Trj.DelAction(Int)
    Trj.Close()
print "timing:", Int.TimeElapsed
print "\n"
Sys.Measures.Reset()
if WriteTraj:
    exit()

Int.Method.Thermostat = Int.Method.ThermostatNone
Sys.Measures.Reset()
print "MD-NVE"
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




