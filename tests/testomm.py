#usr/bin/env python

### Testing for LAMMPS in SIM suite.
### coded by KS
### Example script running the openMM export with a triatomic:
###   1 harmonic bond, 1 rigid bond, external potential, pair interaction
### In the end, compares energies

import time

import numpy as np

#import os, sys
#dir_path = os.path.dirname(os.path.realpath(__file__))
#newpath = "/".join(dir_path.split("/")[:-3])
#sys.path.insert(0,newpath)
#print("expected sim path: {}".format(newpath))
import relpath
relpath.setrelpath(__file__)
import sim
print("actual sim path: {}".format(sim.__file__)) #make sure running correct version

#sim.export.omm.UseTabulated = True
sim.export.omm.NPairPotentialBins = 1000

NStepsMin = 0

#parameter for number of molecules
NMol = 10

#bondlength
bondl = 1.0
bondF = 10.
bondsig = np.sqrt(0.5/bondF)

#setpoint temperature
TempSet = 1.0

np.random.seed(12345)

#define two atom types; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Color = (1,0,0))
AtomTypeB = sim.chem.AtomType("B", Mass = 1.0, Color = (1,0,0))
AtomTypeC = sim.chem.AtomType("C", Mass = 1.0, Color = (1,0,0))

#define a molecule type; give it a name and list of component atom types
MolTypeA = sim.chem.MolType("MA", [AtomTypeA,AtomTypeB,AtomTypeC])
MolTypeA.Bond(0,1) #harmonic bond
MolTypeA.Bond(1,2,bondl) #constrained bond
NMonInMol = 3
COMPos = np.array([ [0,0,0],[2,0,0],[3,0,0] ])
MolTypeA.COMPos = COMPos #need to set COMPos if want to use CubicLattice function to generate position for rigid bonds

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolTypeA], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "test_omm"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolTypeA.New()

#set the system box length sizes
Sys.BoxL = 10

#--- make a new potential energy term ---
P = sim.potential.LJGaussian(Sys, Cut = 3,
                     Filter = sim.atomselect.Pairs,
                     Epsilon = 0.0, Sigma = 1.0, B = 1.0, Kappa=1.0, Dist0=0.0, 
                     Shift = True, Label = "LJG")
P1 = sim.potential.LJGaussian(Sys, Cut = 3,
                     Filter = sim.atomselect.PolyFilter( [AtomTypeA, AtomTypeA] ),
                     Epsilon = 0.0, Sigma = 1.0, B = 1.0, Kappa=1.0, Dist0=0.0, 
                     Shift = True, Label = "LJG")
#add this term to the system forcefield
Sys.ForceField.append(P)
Sys.ForceField.append(P1)

#--- PBonding ---
FilterAB = sim.atomselect.PolyFilter([AtomTypeA,AtomTypeB], Bonded=True)
PBondAB = sim.potential.Bond(Sys, Filter = FilterAB, Dist0 = bondl, FConst = bondF)
Sys.ForceField.extend([PBondAB])
print("expecting bond std ~ {}".format(bondsig))


#--- External ---
PExt = sim.potential.ExternalSinusoid(Sys, Filter = sim.atomselect.PolyFilter([AtomTypeA]),UConst=1.0,NPeriods=1,PlaneAxis=0,PlaneLoc=0.0,Label='ExtSin')
Sys.ForceField.append(PExt)

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
Sys.TempSet = TempSet
Int.Method.Thermostat = Int.Method.ThermostatLangevin
Int.Method.LangevinGamma = 1.0
Int.Method.TimeStep = 0.01
sim.system.positions.CubicLattice(Sys, Random = 0.1)
sim.system.velocities.Canonical(Sys, Temp = TempSet)


Prefix = "testomm_"
ret = sim.export.omm.MakeOpenMMTraj(Sys, Verbose = True,
                                       NStepsMin = NStepsMin,
                                       NStepsEquil = 1000,
                                       NStepsProd = 10000,
                                       WriteFreq = 1000,
                                       Prefix = Prefix,
                                       DelTempFiles = False)
Traj, TrajFile = ret

"""
logdata = np.loadtxt(Prefix + "log.txt")
for (i, Pos) in enumerate(Traj):
    Sys.Pos = Pos
    Sys.ForceField.Eval()
    syspe = Sys.PEnergy
    trjpe = logdata[i,1]
    err = np.abs(syspe - trjpe)
    print i, syspe, trjpe, err, err/np.abs(syspe)
"""

def dist(Pos1,Pos2,BoxL):
    r1 = np.mod(Pos1,BoxL) 
    r2 = np.mod(Pos2,BoxL)
    r = r2 - r1

    r = r - BoxL*np.round(r/BoxL)
    return np.sqrt(np.sum(r**2))

PErelerrs = []
flexbondl = []
flexbondsig = []
rigidbondl = []
for (i, Pos) in enumerate(Traj):
    #print(Pos)
    Sys.Pos = Pos
    Sys.ForceField.Eval()
    syspe = Sys.PEnergy
    trjpe = Traj.FrameData["PEnergy"]
    err = np.abs(syspe - trjpe)
    print("#state  simPE  ommPE  error  relative-error")
    print i, syspe, trjpe, err, err/np.abs(syspe)
    PErelerrs.append(err/np.abs(syspe))

    flexiblebonds = np.zeros(NMol)
    for ii,ind in enumerate(np.arange(0,NMol*NMonInMol,NMonInMol)):
        r1 = Pos[ind,:]
        r2 = Pos[ind+1,:]
        flexiblebonds[ii] = dist(r1,r2,Sys.BoxL) 
        
    rigidbonds = np.zeros(NMol)
    for ii,ind in enumerate(np.arange(1,NMol*NMonInMol,NMonInMol)):
        rigidbonds[ii] = dist(Pos[ind,:],Pos[ind+1,:],Sys.BoxL) 


    flexbondl.append(flexiblebonds.mean())
    flexbondsig.append(flexiblebonds.std())
    rigidbondl.append(rigidbonds.mean())
    print("#<flexible bondl>/bondl, std flexible bond/expected std,    <rigid bondl>/bondl")
    print("{}\t{}\t{}".format(flexiblebonds.mean()/bondl, flexiblebonds.std()/bondsig, rigidbonds.mean()/bondl))

PErelerrs = np.array(PErelerrs)
flexbondl = np.array(flexbondl)
flexbondsig = np.array(flexbondsig)
rigidbondl = np.array(rigidbondl)


print("")
print("=== Summary ===")
print("Relative error in PE: {}".format( PErelerrs.mean() ))
print("Flexible bond relative length, relative std: {}, {}".format(flexbondl.mean()/bondl, flexbondsig.mean()/bondsig))
print("Rigid bond relative length: {}".format(rigidbondl.mean()))


