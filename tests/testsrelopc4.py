#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS, modified by Kevin to map to smeared charge model
#   Todo: why did Scott use Beta=1?
#   Todo: get Lammps Working

import os
import numpy as np
import sim

#THIS FILE COARSE-GRAINS SPCE TO A RIGID TWO SITE DIPOLE MODEL

#parameter for number of molecules
NMol = 267

#density  kg/m^3
Rho = 997.
#convert to number density
Rho = Rho * 6.022e23 / 0.0180153 / 1.e30


#setpoint temperature
TempSet = 298.
kB = 1.380658e-23
kT = kB*TempSet
kTkJmol = kT/1000.0*6.022e23
kTkcalmol = kTkJmol/4.184

#use lammps for CG simulations?
UseLammps = False


#CG parameters -- we fix the geometry of the CG molecule, and the smearing widths.
#                 we let Srel determine the repulsion parameter and partial charges
astd = 1 #A, in the convention where `a` is the stdev of the Gaussian
aev  = 1 #A, in the convention where `a` is the stdev of the Gaussian 
rcut = 3.1658*2.5

np.random.seed(12345)


#CHECK IF A TRAJECTORY EXISTS
trjfn = "{}/sampletraj/opc4.lammpstrj".format(sim.__path__[0])

if os.path.isfile(trjfn):
    Trj = sim.traj.Lammps(trjfn)
    #todo: verify NMol
    #todo: figure out BoxL
    L = 20.0
else:
    raise ValueError('could not locate reference opc4 trajectory')
    
#now make the CG system
print "NOW MAKING CG SYSTEM"
#define atom types
AtomTypeO = sim.chem.AtomType("O", Mass = 15.99, Charge = -1.3582, Color = (1,0,0))
AtomTypeH = sim.chem.AtomType("H", Mass = 2*1.01588, Charge = 1.3582, Color = (0,1,0))
#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("W", [AtomTypeO, AtomTypeH], Rigid = True)   
#define the relative center of mass positions for the rigid molecule
COMPos = np.array([[0.00000,  0.00000,  0.00000],
                   [0.00000,  0.53950,  0.00000]], dtype=float)
#add rigid bonds
bl01 = sim.geom.Length(COMPos[0] - COMPos[1])
MolType.Bond(0,1,RLength = bl01)
MolType.COMPos = COMPos   
#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.AtomicUnits)
#make a system that exists in this world
SysName = "testsrelopc4"
Sys = sim.system.System(World, Name = SysName)  
Sys.ForceField.Globals.Charge.Fixed = False
#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()   
#set the system box length sizes
#L = (float(NMol) / Rho)**(1./3)
Sys.BoxL = L
print "Box dimensions: ", Sys.BoxL
print "Box volume: ", np.prod(Sys.BoxL)    

#=== make a new potential energy term ===
#Filter = sim.atomselect.PolyFilter([AtomTypeO, AtomTypeO])
#P1 = sim.potential.PairSpline(Sys, Cut = 3.1658 * 2.5, NKnot = 25,
#                              Filter = Filter, Label = "SP")
evParam={"cut":rcut, "Sigma":1.0, "B":5.0*kTkcalmol, "Kappa":0.25/aev/aev, "Dist0":0.0}
#   note that sim's convention is LAMMPS convention, B*exp(-kappa*(r/sigma)^2)
#   can either set Sigma=aev, kappa=0.25, or Sigma=1.0, Kappa=0.25/aev^2
P1 = sim.potential.LJGaussian(Sys, Filter=sim.atomselect.Pairs, Cut=evParam["cut"], 
                            Sigma=evParam["Sigma"], Epsilon=0.0, 
                            B=evParam["B"], Kappa=evParam["Kappa"], 
                            Dist0=evParam["Dist0"], Label="ev")
P2 = sim.potential.Ewald(Sys, Label = "EW", ExcludeBondOrd = 3, 
                         Cut = rcut, Shift = True)
P3 = sim.potential.SmearedCoulombEwCorr(Sys, Filter=sim.atomselect.Pairs, Label = "EW-smear-corr", 
                         Cut = rcut, Shift = True, 
                         FixedCoef = True, BornA = astd*np.sqrt(np.pi))


#add this term to the system forcefield
Sys.ForceField.extend([P1, P2, P3])    
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

#make the mapping function for opc4 -- right now place chgs on Oxygen and midpt b/t H. Ignore virtual site.
Map = sim.atommap.PosMap()
for i in range(NMol):
    Map += [sim.atommap.AtomMap(Atoms1 = 4*i, Atom2 = 2*i)]
    Map += [sim.atommap.AtomMap(Atoms1 = [4*i+1, 4*i+2], Atom2 = 2*i+1)]

Opt = sim.srel.OptimizeTrajClass(Sys, Map, Beta = 1., Traj = Trj, UseTarHists = False, 
                                 FilePrefix = "testsrelopc4")
#make sure the overall system charge is zero
Opt.ConstrainNeutralCharge()
#fix most of the LJGaussian Parameters
Sys.ForceField[0].FreezeSpecificParam([0,1,3,4])


if UseLammps:
    Opt = sim.srel.UseLammps(Opt)

#setup the run
Opt.StepsEquil = 5000
Opt.StepsProd = 100000
Opt.StepsStride = 100
Int.Method = Int.Methods.VVIntegrate        
Int.Method.Thermostat = Int.Method.ThermostatNoseHoover
Opt.Run()
