#usr/bin/env python

### Testing for potentials in SIM suite.
### coded by MSS

import sys
import numpy as np
import sim
from sim import chem
from sim import atomselect
from sim import potential
#import sim.system.visualize

np.random.seed(429342)

PosDelta = 0.000001
BoxLDelta = 0.000001
ParamDelta = 0.00000001
ParamMult = 1.00000001
ArgHistNBin = 1000000
ShowDifference1 = False
AbsErrTol = 0.001
FracErrTol = 0.001

ErrList = []

def UpdateErrors(PotLabel, Label, Val1, Val2):
    global ErrList
    AbsErr = abs(Val1-Val2)
    FracErr = AbsErr * 2.0 / abs(Val1 + Val2 + 1.e-300)
    if AbsErr > AbsErrTol and FracErr > FracErrTol:
        ErrList.append((PotLabel, Label, AbsErr, FracErr))



def Test(PClass):
    global ErrList
    
    Name = PClass.__name__
    print "="*70
    print "RUNNING TESTS FOR POTENTIAL %s." % Name

    NAtomPerMolAB = 18
    NMolAB = 9
    NAtomPerMolC = 1
    NMolC = 54
    NMol = NMolAB + NMolC
    MoveAtom = 1
    np.random.seed(12345)
       
    #make a dummy system
    AtomTypeA = chem.AtomType("A", Mass = 1.2, Charge = -0.6)
    AtomTypeB = chem.AtomType("B", Mass = 0.8, Charge = 0.6)
    AtomTypeC = chem.AtomType("C", Mass = 0.6, Charge = 0.0)
    MolTypeAB = chem.MolType("MAB", [AtomTypeA, AtomTypeB]*(NAtomPerMolAB/2))
    for i in range(len(MolTypeAB)-1):
        MolTypeAB.Bond(i, i+1)
    MolTypeC = chem.MolType("MC", [AtomTypeC]*(NAtomPerMolC))
    World = chem.World([MolTypeAB, MolTypeC], Dim = 3, Units = sim.units.DimensionlessUnits)

    SysName = "test" + Name
    Sys = sim.system.System(World, Name = SysName)
    for i in range(NMolAB):
        Sys += MolTypeAB.New()
    for i in range(NMolC):
        Sys += MolTypeC.New()
    
    Sys.BoxL = 5

    #cubic lattice positions
    sim.system.positions.CubicLattice(Sys, Random = 0.15)
    
    if PClass.Type == potential.ptypes.FieldPotential:
        Filter = atomselect.All
    elif PClass.Type == potential.ptypes.PairPotential:
        Filter = atomselect.Pairs
    elif PClass.Type == potential.ptypes.AnglePotential:
        Filter = atomselect.BondTriples
    elif PClass.Type == potential.ptypes.TorsionPotential:
        Filter = atomselect.BondQuartets
    elif PClass.Type == potential.ptypes.LocalDensityPotential:
        Filter = atomselect.Pairs.copy()
        Filter.Ordered = True
    elif PClass.Type == potential.ptypes.GlobalPotential:
        Filter = atomselect.All
    else:
        raise TypeError("Did not recognize type of potential %s" % Name)  

    Potential = PClass(Sys, Label = "Test", Filter = Filter, **PClass.TestArgs)
    
    Potential.Arg.SetupHist(NBin = ArgHistNBin, ReportNBin = 1000)
   
    Potential.SetupTest()
    
    Sys.ForceField.append(Potential)
    if Potential.UsesCharge:
        Sys.ForceField.Globals.Charge.Fixed = False
    Sys.Load()    
    
    
    #check arguments
    print "Running argument test."
    Sys.ForceField.EvalArgStatsThis()
    print "Argument mins  : ", Potential.Arg.Min
    print "Argument maxs  : ", Potential.Arg.Max
    print "Argument avgs  : ", Potential.Arg.Avg
    print "Argument stds  : ", Potential.Arg.Std
    print "Argument count : ", Potential.Arg.Count
    print ""
    
    #system energies
    Sys.ForceField.Eval()
    print "System potential energy =", Sys.PEnergy, "\n"

    #check force
    print "Running force calculation test."
    Atom = Sys.Atom[MoveAtom]
    Sys.Flags.CalcForce = True
    Sys.Flags.CalcVirial = True
    
    Sys.ForceField.Eval()
    PEnergy1 = Sys.PEnergy
    Force1 = Atom.Force
    DiffForce = Force1 * 0.
    for i in range(3):
        Atom.Pos[i] += PosDelta
        Sys.ForceField.Eval()
        PEnergy2 = Sys.PEnergy
        DiffForce[i] = -(PEnergy2 - PEnergy1) / PosDelta
        Atom.Pos[i] -= PosDelta
    Sys.ForceField.Eval()
    Force2 = Atom.Force
    CodeForce = 0.5 * (Force1 + Force2)
    Error = np.sqrt(np.sum((CodeForce - DiffForce)**2))
    print "  Computed force  : %11.5f, %11.5f, %11.5f" % tuple(CodeForce)
    print "  Difference force: %11.5f, %11.5f, %11.5f" % tuple(DiffForce)
    print "  Error           : %11.5f" % Error
    print "\n"
    UpdateErrors(Name, "Force_x", CodeForce[0], DiffForce[0])
    UpdateErrors(Name, "Force_y", CodeForce[1], DiffForce[1])
    UpdateErrors(Name, "Force_z", CodeForce[2], DiffForce[2])
    del Atom
    
    
    #check virial
    print "Running virial test."
    Sys.Flags.CalcVirial = True
    Sys.ForceField.Eval()
    Virial1 = Sys.Virial
    OldBoxL = Sys.BoxL.copy()
    OldPos = Sys.Pos.copy()
    OldVol = np.prod(Sys.BoxL)
    Scale = 1. + BoxLDelta
    Sys.BoxL = OldBoxL * Scale
    Sys.Pos = OldPos * Scale
    Sys.ForceField.Eval()
    PEnergy2a = Sys.PEnergy
    Vol2a = np.prod(Sys.BoxL)
    Scale = 1. - BoxLDelta
    Sys.BoxL = OldBoxL * Scale
    Sys.Pos = OldPos * Scale
    Sys.ForceField.Eval()
    PEnergy2b = Sys.PEnergy
    Vol2b = np.prod(Sys.BoxL)
    Virial2 = 3. * OldVol * (PEnergy2a - PEnergy2b) / (Vol2a - Vol2b)
    Error = np.abs(Virial2 - Virial1)
    print "  Computed virial  : %11.5f" % Virial1
    print "  Difference virial: %11.5f" % Virial2
    print "  Error            : %11.5f" % Error
    print "\n"
    UpdateErrors(Name, "Virial", Virial1, Virial2)
    Sys.BoxL = OldBoxL
    Sys.Pos = OldPos
    Sys.ForceField.Eval()
    
    
    
    #check derivatives wrt parameters
    FF = Sys.ForceField
    for Stage in [1,2]:
        
        if Stage == 1:    
            print "Running parameter derivative calculation test -- explicit calcs."
            Sys.Flags.CalcDUParam = True
            Sys.Flags.CalcDWParam = True
            FF.Eval()
            RefEnergy = Sys.PEnergy
            def FFEval():
                return FF.Eval()
        elif Stage == 2:
            if Potential.NoHistEval: continue
            print "Running parameter derivative calculation test -- histogram calcs."
            Sys.Flags.CalcDUParam = True
            Sys.Flags.CalcDWParam = False
            FF.Arg.ResetStats()
            FF.EvalArgStatsAdd()
            FF.Arg.NormalizeStats()
            Potential.Arg.SetupHist()
            FF.Arg.ResetHist()
            FF.EvalArgHistAdd()
            FF.Arg.NormalizeHist()
            FF.EvalArg()
            HistEnergy = Sys.PEnergy
            print "Explicit energy was: %11.5f" % RefEnergy
            print "Histogram energy is: %11.5f" % HistEnergy
            print "Error: %11.5f" % (abs(RefEnergy-HistEnergy))
            UpdateErrors(Name, "Histogram energy", RefEnergy, HistEnergy)
            file("argumenthistogram.txt","w").write(Potential.Arg.HistString())
            def FFEval():
                return FF.EvalArg()

        for (i, iVal) in enumerate(FF.Param.Val):
            if FF.Param.Fixed[i]: continue
            if Potential.UsesFluct and FF.Param.Names[i].startswith("FLKnots"): continue
            iName = FF.Param.Names[i].replace("Test:","")
            print "Parameter %s" % iName
            iVal1 = iVal
            FFEval()
            CodeDU = FF.Param.DU[i]
            CodeDDU = FF.Param.DDU.Matrix[i,:]
            CodeDW = FF.Param.DW[i]
            CodeDDW = FF.Param.DDW.Matrix[i,:]
            PEnergyi1j1 = Sys.PEnergy
            Viriali1j1 = Sys.Virial
            #change value of i
            if iVal == 0.:
                iVal2 = ParamDelta
                iVal3 = ParamDelta * 2
            else:
                iVal2 = iVal * ParamMult
                iVal3 = iVal * (2*ParamMult - 1)
            #get derivative
            FF.Param.Val[i] = iVal2
            FF.Update()
            FFEval()
            PEnergyi2j1 = Sys.PEnergy
            Viriali2j1 = Sys.Virial
            DiffDUj1 = (PEnergyi2j1 - PEnergyi1j1) / (iVal2 - iVal1)
            DiffDWj1 = (Viriali2j1 - Viriali1j1) / (iVal2 - iVal1)
            CodeDU2 = FF.Param.DU[i]
            CodeDW2 = FF.Param.DW[i]
            print "  Energy first derivative:"
            print "    Computed   : %11.5f" % CodeDU
            print "    Difference : %11.5f   [err: %11.5f]" % (DiffDUj1, abs(CodeDU - DiffDUj1))
            UpdateErrors(Name, "%s DU/D(%s)" % (Name, iName), CodeDU, DiffDUj1)
            if Stage == 1:
                print "  Virial first derivative:"
                print "    Computed   : %11.5f" % CodeDW
                print "    Difference : %11.5f   [err: %11.5f]" % (DiffDWj1, abs(CodeDW - DiffDWj1)) 
                UpdateErrors(Name, "%s DW/D(%s)" % (Name, iName), CodeDW, DiffDWj1)
            #loop through cross-derivatives
            for (j, jVal) in enumerate(FF.Param.Val):
                if FF.Param.Fixed[j]: continue
                jName = FF.Param.Names[j].replace("Test:","")
                jVal1 = jVal
                if i==j:
                    FF.Param.Val[i] = iVal3
                    FF.Update()
                    FFEval()
                    PEnergyi3j1 = Sys.PEnergy
                    Viriali3j1 = Sys.Virial
                    DiffDDU = (PEnergyi3j1 - 2*PEnergyi2j1 + PEnergyi1j1) / ((iVal2 - iVal1)*(iVal3- iVal2))
                    DiffDDU2 = (CodeDU2 - CodeDU) / (iVal2 - iVal1)
                    DiffDDW = (Viriali3j1 - 2*Viriali2j1 + Viriali1j1) / ((iVal2 - iVal1)*(iVal3- iVal2))
                    DiffDDW2 = (CodeDW2 - CodeDW) / (iVal2 - iVal1)
                else:
                    #change value of j
                    if jVal == 0.:
                        jVal2 = ParamDelta
                    else:
                        jVal2 = jVal * ParamMult
                    FF.Param.Val[i] = iVal1
                    FF.Param.Val[j] = jVal2
                    FF.Update()
                    FFEval()
                    DiffDDU2 = (FF.Param.DU[i] - CodeDU) / (jVal2 - jVal1)
                    DiffDDW2 = (FF.Param.DW[i] - CodeDW) / (jVal2 - jVal1)
                    PEnergyi1j2 = Sys.PEnergy
                    Viriali1j2 = Sys.Virial
                    FF.Param.Val[i] = iVal2
                    FF.Param.Val[j] = jVal2
                    FF.Update()
                    FFEval()
                    PEnergyi2j2 = Sys.PEnergy
                    Viriali2j2 = Sys.Virial
                    DiffDUj2 = (PEnergyi2j2 - PEnergyi1j2) / (iVal2 - iVal1)
                    DiffDWj2 = (Viriali2j2 - Viriali1j2) / (iVal2 - iVal1)
                    #need an all close HERE?
                    DiffDDU = (DiffDUj2 - DiffDUj1) / (jVal2 - jVal1)
                    DiffDDW = (DiffDWj2 - DiffDWj1) / (jVal2 - jVal1)
                print "  Energy second derivative with %s:" % jName
                print "    Computed   : %11.5f" % CodeDDU[j]
                if ShowDifference1:
                    print "    Difference : %11.5f   [err: %11.5f]" % (DiffDDU, abs(CodeDDU[j] - DiffDDU))
                print "    Difference2: %11.5f   [err: %11.5f]" % (DiffDDU2, abs(CodeDDU[j] - DiffDDU2))
                UpdateErrors(Name, "%s D2U/D(%s)D(%s)" % (Name, iName, jName), CodeDDU[j], DiffDDU2)
                if Stage == 1:
                    print "  Virial second derivative with %s:" % jName
                    print "    Computed   : %11.5f" % CodeDDW[j]
                    if ShowDifference1:
                        print "    Difference : %11.5f   [err: %11.5f]" % (DiffDDW, abs(CodeDDW[j] - DiffDDW))
                    print "    Difference2: %11.5f   [err: %11.5f]" % (DiffDDW2, abs(CodeDDW[j] - DiffDDW2))
                    UpdateErrors(Name, "%s D2W/D(%s)D(%s)" % (Name, iName, jName), CodeDDW[j], DiffDDW2)
                #go back to original
                FF.Param.Val[j] = jVal1
                FF.Param.Val[i] = iVal1
                FF.Update()
            #go back to original
            FF.Param.Val[i] = iVal1
            FF.Update()
            print "\n"
            
    if Potential.UsesFluct:
    
        print "Running parameter derivative calculation test -- fluctuations.\n"
        
        def FFEval():
            return FF.Eval()
        
        #test the fluctuation calculations     
        Sys.Flags.CalcDUParam = False
        Sys.Flags.CalcDWParam = False
        Sys.Flags.CalcFluct = True
        Sys.ForceField.FluctE - Sys.Units.EScale * 2.
        Sys.ForceField.FluctBeta = 0.8 
  
        for (i, iVal) in enumerate(FF.Param.Val):
            if FF.Param.Fixed[i]: continue
            iName = FF.Param.Names[i].replace("Test:","")
            print "Parameter %s" % iName
            iVal1 = iVal
            FFEval()
            CodeDU = FF.Param.DU[i]
            CodeDDU = FF.Param.DDU.Matrix[i,:]
            NLnPi1j1 = FF.FluctTerm
            #change value of i
            if iVal == 0.:
                iVal2 = ParamDelta
                iVal3 = ParamDelta * 2
            else:
                iVal2 = iVal * ParamMult
                iVal3 = iVal * (2*ParamMult - 1)
            #get derivative
            FF.Param.Val[i] = iVal2
            FF.Update()
            FFEval()
            NLnPi2j1 = FF.FluctTerm
            DiffDUj1 = (NLnPi2j1 - NLnPi1j1) / (iVal2 - iVal1)
            CodeDU2 = FF.Param.DU[i]
            print "  -lnP first derivative:"
            print "    Computed   : %11.5f" % CodeDU
            print "    Difference : %11.5f   [err: %11.5f]" % (DiffDUj1, abs(CodeDU - DiffDUj1))
            UpdateErrors(Name, "%s DNlnP/D(%s)" % (Name, iName), CodeDU, DiffDUj1)
            #loop through cross-derivatives
            for (j, jVal) in enumerate(FF.Param.Val):
                if FF.Param.Fixed[j]: continue
                jName = FF.Param.Names[j].replace("Test:","")
                jVal1 = jVal
                if i==j:
                    FF.Param.Val[i] = iVal3
                    FF.Update()
                    FFEval()
                    NLnPi3j1 = FF.FluctTerm
                    DiffDDU = (NLnPi3j1 - 2*NLnPi2j1 + NLnPi1j1) / ((iVal2 - iVal1)*(iVal3- iVal2))
                    DiffDDU2 = (CodeDU2 - CodeDU) / (iVal2 - iVal1)
                else:
                    #change value of j
                    if jVal == 0.:
                        jVal2 = ParamDelta
                    else:
                        jVal2 = jVal * ParamMult
                    FF.Param.Val[i] = iVal1
                    FF.Param.Val[j] = jVal2
                    FF.Update()
                    FFEval()
                    DiffDDU2 = (FF.Param.DU[i] - CodeDU) / (jVal2 - jVal1)
                    NLnPi1j2 = FF.FluctTerm
                    FF.Param.Val[i] = iVal2
                    FF.Param.Val[j] = jVal2
                    FF.Update()
                    FFEval()
                    NLnPi2j2 = FF.FluctTerm
                    DiffDUj2 = (NLnPi2j2 - NLnPi1j2) / (iVal2 - iVal1)
                    #need an all close HERE?
                    DiffDDU = (DiffDUj2 - DiffDUj1) / (jVal2 - jVal1)
                print "  -lnP second derivative with %s:" % jName
                print "    Computed   : %11.5f" % CodeDDU[j]
                if ShowDifference1:
                    print "    Difference : %11.5f   [err: %11.5f]" % (DiffDDU, abs(CodeDDU[j] - DiffDDU))
                print "    Difference2: %11.5f   [err: %11.5f]" % (DiffDDU2, abs(CodeDDU[j] - DiffDDU2))
                UpdateErrors(Name, "%s D2NlnP/D(%s)D(%s)" % (Name, iName, jName), CodeDDU[j], DiffDDU2)
                #go back to original
                FF.Param.Val[j] = jVal1
                FF.Param.Val[i] = iVal1
                FF.Update()
            #go back to original
            FF.Param.Val[i] = iVal1
            FF.Update()
            print "\n"        
            
    #atom monte carlo routines
    Sys.Flags.CalcsOff()
    print ""
    
    #deletion of a molecule
    print "Running molecule deletion test..."
    Sys.ForceField.Eval()
    Sys.ForceField.SaveEneState()
    PEnergy0 = Sys.PEnergy        
    Sys.TargetMol = 0
    Sys.ForceField.Eval(Mode = -2)
    PEnergy1 = Sys.PEnergy
    Sys.MolActive[0] = -1        
    Sys.ForceField.Eval()        
    PEnergy2 = Sys.PEnergy
    Error = np.abs(PEnergy1 - PEnergy2)
    print "  Molecule update energy: %11.5f" % PEnergy1
    print "  Full computed energy  : %11.5f" % PEnergy2
    print "  Error                 : %11.5f" % Error
    UpdateErrors(Name, "Mol deletion ene", PEnergy2, PEnergy1)
    Sys.MolActive[0] = 1
    Sys.ForceField.RevertEneState()
    PEnergy3 = Sys.PEnergy
    Sys.ForceField.Eval()
    PEnergy4 = Sys.PEnergy
    Error = np.abs(PEnergy4 - PEnergy3)
    print "  Energy after reverting: %11.5f" % PEnergy4
    print "  Full computed energy  : %11.5f" % PEnergy3
    print "  Error                : %11.5f" % Error        
    UpdateErrors(Name, "Mol del/revert ene", PEnergy4, PEnergy3)
    print ""

    #insertion of a molecule
    print "Running molecule insertion test..."
    Sys.MolActive[0] = -1
    Sys.ForceField.Eval()
    Sys.ForceField.SaveEneState()
    PEnergy0 = Sys.PEnergy   
    Sys.TargetMol = 0        
    Sys.ForceField.Eval(Mode = +2)
    PEnergy1 = Sys.PEnergy
    Sys.MolActive[0] = 1        
    Sys.ForceField.Eval()        
    PEnergy2 = Sys.PEnergy
    Error = np.abs(PEnergy1 - PEnergy2)
    print "  Molecule update energy: %11.5f" % PEnergy1
    print "  Full computed energy  : %11.5f" % PEnergy2
    print "  Error                 : %11.5f" % Error
    UpdateErrors(Name, "Mol insertion ene", PEnergy2, PEnergy1)
    Sys.MolActive[0] = -1
    Sys.ForceField.RevertEneState()
    PEnergy3 = Sys.PEnergy
    Sys.ForceField.Eval()
    PEnergy4 = Sys.PEnergy
    Error = np.abs(PEnergy4 - PEnergy3)
    print "  Energy after reverting: %11.5f" % PEnergy4
    print "  Full computed energy  : %11.5f" % PEnergy3
    print "  Error                 : %11.5f" % Error        
    UpdateErrors(Name, "Mol ins/revert ene", PEnergy4, PEnergy3)  
    Sys.MolActive[0] = 1
    Sys.ForceField.Eval()
    print ""
            
    #displacement of a molecule
    print "Running molecule displacement test..."
    for i in [0, 1]:
        Sys.ForceField.Eval()
        Sys.ForceField.SaveEneState()
        PEnergy0 = Sys.PEnergy   
        Sys.TargetAtom = i        
        Sys.ForceField.Eval(Mode = -1)
        OldiPos = Sys.Pos[i].copy()
        Sys.Pos[i] = Sys.Pos[i] + 0.1 * (np.random.rand(3) - 0.5)    
        Sys.ForceField.Eval(Mode = +1) 
        PEnergy1 = Sys.PEnergy
        Sys.ForceField.Eval()
        PEnergy2 = Sys.PEnergy
        Error = np.abs(PEnergy1 - PEnergy2)
        print "  Molecule update energy: %11.5f" % PEnergy1
        print "  Full computed energy  : %11.5f" % PEnergy2
        print "  Error                 : %11.5f" % Error
        UpdateErrors(Name, "Atom displace ene", PEnergy2, PEnergy1)
        Sys.Pos[i] = OldiPos
        Sys.ForceField.RevertEneState()
        PEnergy3 = Sys.PEnergy
        Sys.ForceField.Eval()
        PEnergy4 = Sys.PEnergy
        Error = np.abs(PEnergy4 - PEnergy3)
        print "  Energy after reverting: %11.5f" % PEnergy4
        print "  Full computed energy  : %11.5f" % PEnergy3
        print "  Error                 : %11.5f" % Error        
        UpdateErrors(Name, "Atom displace/revert ene", PEnergy4, PEnergy3)  
        Sys.ForceField.Eval()
    print ""
               
    Sys.Cleanup()
    World.Cleanup()
    

#make a list of all potential classes
def GetAllSubclasses(x):
  l = list(x.__subclasses__())
  for y in l:
    l.extend(GetAllSubclasses(y))
  l = [y for (i,y) in enumerate(l) if not y in l[i+1:]]
  return l


PClasses = GetAllSubclasses(potential.PotentialClass)


if len(sys.argv) > 1:
    Names = [x.__name__ for x in PClasses]
    for Arg in sys.argv[1:]:
        if not Arg in Names:
            print "Could not find class definition for %s" % Arg
    PClasses = [x for x in PClasses if x.__name__ in sys.argv[1:]]

    
for PClass in PClasses:
    if not hasattr(PClass, "TestArgs"):
        print "Could not find testing arguments for class:", str(PClass)
        continue
    Test(PClass)
            
print "\n" + "="*70
print "SUMMARY OF ERRORS"
print "  %-50s %-12s %-12s" % ("value", "abs err", "frac err")
if len(ErrList) == 0:
    print "-----NONE-----"
OldPLabel = None
for (PLabel, Label, AbsErr, FracErr) in ErrList:
    if OldPLabel != PLabel:
        print PLabel
        OldPLabel = PLabel
    if len(Label) > 50:
        print "  %s" % Label
        print "  %-50s %12.5e %12.5e" % ("", AbsErr, FracErr)   
    else:
        print "  %-50s %12.5e %12.5e" % (Label, AbsErr, FracErr)              
            
            
            
    
    
    