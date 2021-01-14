#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS

import chem
import units
import system
import potential
import atomselect
import srel
import traj
import numpy as np

np.random.seed(3423)

AtomTypeH = chem.AtomType("H", Charge = 0.4238)
AtomTypeO = chem.AtomType("O", Charge = -0.8476)
MolType = chem.MolType("M", [AtomTypeO, AtomTypeH, AtomTypeH])
World = chem.World([MolType], Dim = 3, Units = units.AtomicUnits)

SysName = "srelw"  
Sys = system.System(World, Name = SysName)

for i in range(216):
    Sys += MolType.New()
Sys.BoxL[:] = (216 * 0.01801528 / 1000 / 6.02e23)**(1./3.) / 1.e-10

OOFilter = (atomselect.Filter(Types = AtomTypeO)
            + atomselect.Filter(Types = AtomTypeO))

if 1 == 1:
    P1 = potential.LJ(Sys, Label = "LJ",
                      Cut = 7.915,
                      Filter = OOFilter,
                      Epsilon = 0.1554,
                      Sigma = 3.166,
                      Fixed = False, 
                      Shift = True)
else:
    P1 = potential.PairSpline(Sys, Label = "PS",
                              Cut = 7., NKnot = 20,
                              Filter = OOFilter)
    P1.EmulateLJ(Sigma=2.7, Epsilon=5)
    
Sys.ForceField.append(P1)

if 1==2:
    Sys.ForceField.Globals.Charge.Fixed = False
    Sys.ForceField.ChargePerAtom = 0.
    ChargeFilter = atomselect.Pairs.copy()
    ChargeFilter.Intra = False
    P2 = potential.ScreenedCoulomb(Sys, Label = "SC",
                                   Cut = 9.0,
                                   Filter = ChargeFilter,
                                   Fixed = False,
                                   FixedKappa = True,
                                   Kappa = 0.1,
                                   Shift = True)
    Sys.ForceField.append(P2)

Sys.Lib.ForceCompile = True
Sys.Load()

Map = atommap.PosMap()
for (i, a) in enumerate(Sys.Atom):
    Map += [atommap.AtomMap(Atoms1 = i, Atom2 = a)]


Trj = traj.Simple("z:/res/srel/spcetraj/r1000T260.trj.bz2")[:1000]
Trj.DataScale["PEnergy"] = 6.947695e-21

Beta = 1. / (260.* Sys.Units.kB)

#make an integrator
Int = integrate.Integrator(Sys)
#choose a method
Int.AddMethod(integrate.VVIntegrate, TimeStep = 0.001)

Optimizer = srel.OptimizeTrajClass(Sys, Map, Beta, Trj)
Optimizer.OutFile = "srelw.txt"
Optimizer.Optimize()


for x in np.arange(2.0,10.0,0.05):
    print x, P1.Val(x)



    