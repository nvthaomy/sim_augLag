#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS

import os
import numpy as np
import sim

#THIS FILE COARSE-GRAINS SPCE TO A RIGID TWO SITE DIPOLE MODEL

#parameter for number of molecules
NMol = 150

#density  kg/m^3
Rho = 1000.
#convert to number density
Rho = Rho * 6.02e23 / 0.018 / 1.e30

#setpoint temperature
TempSet = 298.

#use lammps for CG simulations?
UseLammps = False

np.random.seed(12345)


#CHECK IF A TRAJECTORY EXISTS
trjfn = "{}/sampletraj/spce1.lammpstraj.bz2".format(sim.__path__[0])

if os.path.isfile(trjfn):
    
    Trj = sim.traj.Lammps(trjfn)
    
else:
    
    #make a system and then make the trajectory
    print "FIRST MAKING REFERENCE SPCE TRAJECTORY"
    #define atom types
    AtomTypeO = sim.chem.AtomType("O", Mass = 15.99, Charge = -0.8476, Color = (1,0,0))
    AtomTypeH = sim.chem.AtomType("H", Mass = 1.00794, Charge = 0.4238, Color = (0,1,0))
    #define a molecule type; give it a name and list of component atom types
    MolType = sim.chem.MolType("W", [AtomTypeO, AtomTypeH, AtomTypeH], Rigid = True)   
    #define the relative center of mass positions for the rigid molecule
    COMPos = np.array([[0.00000,  0.00000,  0.00000],
                       [0.81649,  0.57736,  0.00000], 
                       [-0.81649,  0.57736,  0.00000]], dtype=float)
    #add rigid bonds
    bl01 = sim.geom.Length(COMPos[0] - COMPos[1])
    bl02 = sim.geom.Length(COMPos[0] - COMPos[2])
    bl12 = sim.geom.Length(COMPos[1] - COMPos[2])
    MolType.Bond(0,1,RLength = bl01)
    MolType.Bond(0,2,RLength = bl02)
    MolType.Bond(1,2,RLength = bl12)
    MolType.COMPos = COMPos   
    #define the world in terms of a list of molecular species and dimensionality 
    World = sim.chem.World([MolType], Dim = 3, Units = sim.units.AtomicUnits)
    #make a system that exists in this world
    SysName = "testsrelspce1"
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
    #compile and load the system
    Sys.Load()   
    #set initial positions and velocities
    sim.system.positions.CubicLattice(Sys)
    sim.system.velocities.Canonical(Sys, Temp = TempSet)    
    Sys.TempSet = TempSet
    print "energy minimizing..."
    #set up the integrator 
    Int = Sys.Int
    Int.Method = Int.Methods.VVQuench
    Int.Run(1000)        
    #add verbose output of property averages
    Sys.Measures.VerboseOutput(StepFreq = 1000)
    #make a trajectory to recort
    Int.Method = Sys.Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
    print "running MD equilibration..."
    Int.Run(20000) 
    print "timing:", Int.TimeElapsed
    Trj = sim.traj.LammpsWrite(trjfn)
    Trj.AddAction(Sys.Int, StepFreq = 1000)
    print "running MD production..."
    Int.Run(1000000)
    print "timing:", Int.TimeElapsed
    Trj.Close()

#now make the CG system
print "NOW MAKING CG SYSTEM"
#define atom types
AtomTypeO = sim.chem.AtomType("O", Mass = 15.99, Charge = -0.8476, Color = (1,0,0))
AtomTypeH = sim.chem.AtomType("H", Mass = 1.01588, Charge = 0.8476, Color = (0,1,0))
#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("W", [AtomTypeO, AtomTypeH], Rigid = True)   
#define the relative center of mass positions for the rigid molecule
COMPos = np.array([[0.00000,  0.00000,  0.00000],
                   [0.00000,  0.57736,  0.00000]], dtype=float)
#add rigid bonds
bl01 = sim.geom.Length(COMPos[0] - COMPos[1])
MolType.Bond(0,1,RLength = bl01)
MolType.COMPos = COMPos   
#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.AtomicUnits)
#make a system that exists in this world
SysName = "testsrelspce2"
Sys = sim.system.System(World, Name = SysName)  
Sys.ForceField.Globals.Charge.Fixed = False
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
P1 = sim.potential.PairSpline(Sys, Cut = 3.1658 * 2.5, NKnot = 25,
                              Filter = Filter, Label = "SP")
P2 = sim.potential.Ewald(Sys, Label = "EW", ExcludeBondOrd = 3, 
                         Cut = 3.1658 * 2.5, Shift = True)
#add this term to the system forcefield
Sys.ForceField.extend([P1, P2])    
#set up the histograms
for P in Sys.ForceField:
    P.Arg.SetupHist(NBin = 10000, ReportNBin=100)
#compile and load the system
Sys.Load()   
Int = Sys.Int

#set initial positions and velocities
sim.system.positions.CubicLattice(Sys)
sim.system.velocities.Canonical(Sys, Temp = TempSet)    
Sys.TempSet = TempSet

#make the mapping function
Map = sim.atommap.PosMap()
for i in range(NMol):
    Map += [sim.atommap.AtomMap(Atoms1 = 3*i, Atom2 = 2*i)]
    Map += [sim.atommap.AtomMap(Atoms1 = [3*i+1, 3*i+2], Atom2 = 2*i+1)]
    
Opt = sim.srel.OptimizeTrajClass(Sys, Map, Traj = Trj, UseTarHists = False, 
                                 FilePrefix = "testsrelspce")
#make sure the overall system charge is zero
Opt.ConstrainNeutralCharge()

if UseLammps:
    Opt = sim.srel.UseLammps(Opt)

#setup the run
Opt.StepsEquil = 5000
Opt.StepsProd = 100000
Opt.StepsStride = 100
Int.Method = Int.Methods.VVIntegrate        
Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
Opt.Run()
