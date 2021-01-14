#usr/bin/env python

### Testing external sinusoid potential in SIM suite.
### coded by KS

import time

import numpy as np

import sim 
print("sim version: {}".format(sim.__file__))

np.random.seed(12345)
NMol = 1

#define two atom types; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Charge = +1.0, Color = (1,0,0))
AtomTypeB = sim.chem.AtomType("B", Mass = 1.0, Charge = -1.0, Color = (0,1,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA, AtomTypeB])

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "test_external_sinusoid"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()

#set the system box length sizes
Sys.BoxL = [10.0, 10.0, 10.0]

#make a new potential energy term
Filter = sim.atomselect.PolyFilter([AtomTypeA])
Ext = {"UConst":10.0, "NPeriods":1.0, "PlaneAxis":1, "PlaneLoc": 0.0}
P = sim.potential.ExternalSinusoid(Sys, Filter=Filter, UConst=Ext["UConst"], NPeriods=Ext["NPeriods"], PlaneAxis=Ext["PlaneAxis"], PlaneLoc=Ext["PlaneLoc"])
#P = sim.potential.GravityField(Sys,Filter=Filter,GravConst=10)     

#add this term to the system forcefield
Sys.ForceField.append(P)
#NOTE: all potential energy terms must be added before compilation

#compile and load the system
Sys.Load()

#set initial positions and velocities
def getEnergy(x1,x2):
    """returns energy, x1 and x2 should both be 3d lists [a,b,c]"""
    rigidUnits = [x for x in Sys.RigidUnitsIter()] #Iterate over c.o.m., see system/__init__.py
    Sys.SetPos(rigidUnits[0],x1)
    Sys.SetPos(rigidUnits[1],x2)

    Sys.Flags.CalcForce = True
    Sys.Flags.CalcVirial = True
    Sys.Flags.CalcDUParam = True
    Sys.Flags.CalcDWParam = True
    Sys.Flags.CalcForce = True
    Sys.ForceField.Eval()
    Sys.Flags.CalcsOff()
    Sys.Flags.CalcVirial = False
    PEnergy1 = Sys.PEnergy
    Virial1 = Sys.Virial

    r = x1[Ext["PlaneAxis"]]
    print("Energy (Python): {}".format(Sys.ForceField[0].Val(r)) )
    print("Energy (Fortran):{}".format(Sys.PEnergy))

    print("Force (Python): {}".format(-Sys.ForceField[0].DVal(r)))
    print("Force (Fortran):{}".format(Sys.Force[0]))

    tol = 1e-6
    return np.abs( Sys.ForceField[0].Val(r) - Sys.PEnergy ) < tol, Sys.PEnergy, Sys.ForceField[0].Val(r)



tol = 1e-6
testcounter = [0,0]
def countTest(condition,testcounter):
    """testcounter = [#passes, #tests]. In the future should make this an object."""
    if condition:
        testcounter[0] += 1
    else:
        print("failed a test")
    testcounter[1] += 1 



print("\n=== Testing second particle type doesn't change energy, because it wasn't in Filter ===")
x1,x2 = [0,0,0], [3.2,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( np.abs(EFortran) < tol and np.abs(EPython) < tol, testcounter )

x1,x2 = [0,0,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( np.abs(EFortran) < tol and np.abs(EPython) < tol, testcounter )


print("\n=== Testing potential axis: only y-axis should matter ===")
x1,x2 = [0,0,1], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( np.abs(EFortran) < tol and np.abs(EPython) < tol, testcounter )

x1,x2 = [0,1,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( np.abs(EFortran) > tol and np.abs(EPython) > tol, testcounter )

x1,x2 = [1,0,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( np.abs(EFortran) < tol and np.abs(EPython) < tol, testcounter )


print("\n=== Testing sinusoid: make sure fortran matches python ===")
x1,x2 = [0,0,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( match, testcounter )

x1,x2 = [0,1.25,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( match, testcounter )

x1,x2 = [0,2.5,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( match, testcounter )

x1,x2 = [0,3.75,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( match, testcounter )

x1,x2 = [0,5,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( match, testcounter )

x1,x2 = [0,7.5,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( match, testcounter )

x1,x2 = [0,8.75,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( match, testcounter )

x1,x2 = [0,10.0,0], [5.3,5.7,8.9]
print("--> coords: {},{} <--".format(x1,x2))
match,EFortran,EPython = getEnergy(x1, x2)
countTest( match, testcounter )


print("\n=== Summary ===")
print("passed {}/{} tests".format(testcounter[0],testcounter[1]))


'''
print "-", P.RPEnergy, P.SPEnergy, P.RVirial, Sys.PEnergy
Sys.TargetMol = 0
Sys.MolActive[0] = -1
Sys.ForceField.Eval(Mode = -2)
PEnergy2 = Sys.PEnergy
Virial2 = Sys.Virial
print "-", P.RPEnergy, P.SPEnergy, P.RVirial, Sys.PEnergy
Sys.ForceField.Eval()
PEnergy3 = Sys.PEnergy
Virial3 = Sys.Virial
print "-", P.RPEnergy, P.SPEnergy, P.RVirial, Sys.PEnergy
print PEnergy1, PEnergy2, PEnergy3
print Virial1, Virial2, Virial3
'''

