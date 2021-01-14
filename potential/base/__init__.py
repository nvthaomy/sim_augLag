#/usr/bin/env python

### Class definitions for force fields in SIM suite.
### coded by MSS

import numpy as np

import sim.fortran as fortran

import potentialarray as PA
import argarray as AA
import basefortran
import argfortran
import potentialtypes as ptypes

from potentialclass import PotentialClass
from splinepotentialclass import SplinePotentialClass
from splinefluctpotentialclass import SplineFluctPotentialClass
           
        

class ForceField(list):
    """Class for collections of potential energy functions."""

    ModuleVars = ["Terms", "Cut", "CutSq", "OldPEnergy", "OldTerms",
                  "FluctE", "FluctE0", "FluctA0", "FluctBeta", "FluctTerm"]
    ModuleConsts = ["NTerm"]
   
    def __init__(self, Sys, Potentials = None,
                 FixedCharge = True):
        list.__init__(self)
        list.__setattr__(self, "LibVars", {})
        self.Sys = Sys
        if Potentials is None: Potentials = []
        self.extend(Potentials)
        PreCallCode = "Sys = self.Sys"
        #fortran subroutine for energy
        self.EneSub = fortran.Subroutine("calcenergyforces", 
                                         SubStr = "self.EneSub", PreCallCode = PreCallCode)
        self.EneSub.PreLoad = self.PreLoad
        self.EneSub.PostLoad = self.PostLoad
        self.Sys.Lib.append(self.EneSub)
        self.SaveEneSub = fortran.Subroutine("saveenergystate",
                                             SubStr = "self.SaveEneSub", PreCallCode = PreCallCode)
        self.Sys.Lib.append(self.SaveEneSub)     
        self.RevertEneSub = fortran.Subroutine("revertenergystate",
                                               SubStr = "self.RevertEneSub", PreCallCode = PreCallCode)
        self.Sys.Lib.append(self.RevertEneSub)            
        #fortran subroutine for argument distributions
        self.ArgStatsSub = fortran.Subroutine("calcargstats", 
                                              SubStr = "self.ArgStatsSub", PreCallCode = PreCallCode)
        self.Sys.Lib.append(self.ArgStatsSub)
        self.ArgHistSub = fortran.Subroutine("calcarghist",
                                             SubStr = "self.ArgHistSub", PreCallCode = PreCallCode)
        self.Sys.Lib.append(self.ArgHistSub)
        self.ArgEvalSub = fortran.Subroutine("calcargeval",
                                             SubStr = "self.ArgEvalSub", PreCallCode = PreCallCode)
        self.Sys.Lib.append(self.ArgEvalSub)         
        #make global parameters
        self.Globals = PA.Parray("Global")
        Charge = [a.Charge for a in self.Sys.World.AtomTypes]
        self.Globals.Add("Charge", len(Charge), Value = Charge,
                         Fixed = FixedCharge) 
        
        #parameters for energy fluctuations in CG models
        self.FluctE = 0.
        self.FluctE0 = 0.
        self.FluctA0 = 0.
        self.FluctBeta = 0.
        self.FluctTerm = 0.
        
        
    def __str__(self):
        return "ForceField[" + ",".join([str(x) for x in self]) + "]"
 
    def __setattr__(self, name, val):
        Var = self.LibVars.get(name, None)
        if Var is None:
            v = self.__dict__.get(name, None)
            if isinstance(v, np.ndarray):
                v[:] = val
            else:
                object.__setattr__(self, name, val)   
        else:
            if Var.shape:
                Var[:] = val
            else:
                Var.itemset(val)             
                    
    def __getattr__(self, name):
        Var = self.LibVars.get(name, None)
        if Var is None:
            raise AttributeError(name) 
        else:
            if not Var.shape: Var = Var.item()
            return Var
            
    def Cleanup(self):
        for P in self:
            P.Cleanup()
        self.Sys = None
        self.Param = None

    def ArraySizes(self):
        """Returns the sizes of the Param/DUParam and DDUParam arrays."""
        return PA.MasterSizes(self.Globals, self)

    def __CollectArrays(self):
        "Collects all parameters into a master list."
        #check for non-overlap with global names
        for P in self:
            for Name in P.Maps.keys():
                if Name in self.Globals.Maps:
                    raise NameError("%s in %s is already a global." % (Name, P.Name))
                
        #link to masters
        ParrayList = [P.Param for P in self]
        self.Param = PA.MasterParray(self.Globals, ParrayList,
                                     Name = "ForceField")

        #collect all cutoffs
        self.Cut = []
        self.CutSq = []
        for P in self:
            if P.Cut is None:
                self.Cut.append(1.e300)
                self.CutSq.append(1.e300)
            else:
                self.Cut.append(P.Cut)
                self.CutSq.append(P.Cut * P.Cut)
        self.Cut = np.array(self.Cut, float)
        self.CutSq = np.array(self.CutSq, float)

        #collect args
        self.Arg = AA.MasterArgArray([P.Arg for P in self])

        #make a terms array
        self.NTerm = len(self)
        self.Terms = np.zeros(len(self), float)
        self.OldTerms = np.zeros_like(self.Terms)
        self.OldPEnergy = 0.

        #flag for histogram eval
        self.NoHistEval = any([P.NoHistEval for P in self])

    
    def GetUniquePotentials(self, AllTypes = True):
        """Returns a list of (Potential, PotentialName, TypeInd) for all unique potentials."""
        Potentials = []
        for P in self.Sys.ForceField:
            if AllTypes:
                typeindmax = P.Arg.NType
            else:
                typeindmax = 1
            for typeind in range(0,typeindmax):
                if P.Arg.NType == 1:
                    PName = P.Name
                else:
                    PName = "%s(%s)" % (P.Name, P.Arg.TypeLabels[typeind]) 
                Potentials += [(P, PName, typeind)]
        return Potentials
                
                
    def ParamString(self):
        """Returns a string of parameters for this energy function."""
        s = ""
        for P in self:
            s += ">>> POTENTIAL %s\n" % P.Name
            s += P.ParamString()
            s += "\n"
        return s
            
    def SetParamString(self, s, Check = False):
        """Sets parameters from an input string."""
        import ast
        s = s.strip()
        for dat in s.split(">>> POTENTIAL")[1:]:
            ThisName = dat.split("\n")[0].strip()
            ThisData = dat[dat.index("\n")+1:].strip()
            ThisP = None
            for P in self:
                if P.Name == ThisName: 
                    ThisP = P
                    break
            if ThisP is None:
                if Check:
                    raise ValueError("Could not find a potential %s." % ThisName)    
            else:
                P.SetParam(**ast.literal_eval(ThisData))
            
    def SetParam(self, ParamArray, ForceFixed = False):
        """Sets all float parameters for this potential."""
        if self.Param is None:
            raise SystemError("Force field has not been loaded yet.")
        Fixed = self.Param.Fixed
        a = ParamArray.copy()
        #keep fixed values the same
        if not ForceFixed:
            a[Fixed] = self.Param.Val[Fixed]
        #update values
        self.Param.Val = a
        self.Update()

    def Update(self):
        """Updates all potentials for any parameter changes."""
        for P in self:
            P.Update()           

    def PreLoad(self):
        """Run before fortran compilation."""
        #put all parameters into a master array
        self.__CollectArrays()

        #update potentials
        self.Update()
        
        #prepare module-level vars
        self.Sys.Lib.VarPreLoad(self, self.ModuleVars, self.ModuleConsts)
        self.Param.PreLoad(self.Sys.Lib)
        
        #run pre-load command for all potentials
        for (ind, P) in enumerate(self):
            P.FortPrefix = "P%d_" % ind
            P.PreLoad() 
            P.Arg.FortPrefix = "P%d_Arg" % ind
            P.Arg.PreLoad(self.Sys.Lib)

        #main source for energy routine
        fc, fc1, fc2 = basefortran.GetFortCode(self)
        self.EneSub.AddFortCode(fc)
        #source for saving and reverting energies
        self.SaveEneSub.AddFortCode(fc1)
        self.RevertEneSub.AddFortCode(fc2)

        #main source for argument stats routine
        fc = argfortran.GetFortCodeStats(self)
        self.ArgStatsSub.AddFortCode(fc)        

        #main source for argument hist routine
        fc = argfortran.GetFortCodeHist(self)
        self.ArgHistSub.AddFortCode(fc)
        
        #main source for argument hist eval
        fc = argfortran.GetFortCodeEval(self)
        self.ArgEvalSub.AddFortCode(fc)
            

    def PostLoad(self):
        """Run after fortran compilation."""
        #prepare module-level vars
        self.Sys.Lib.VarPostLoad(self)
        self.Param.PostLoad(self.Sys.Lib)
        
        #run post-load command for all potentials
        for P in self:
            P.PostLoad()
            P.Arg.PostLoad(self.Sys.Lib)
       
    def UpdateTerms(self):
        """Updates energy terms."""
        for (i, P) in enumerate(self):
            P.PEnergy = self.Terms[i]

    def Eval(self, Mode = 0):
        """Evaluates system energy and forces."""
        exec(self.EneSub.CallObj)

    def EvalArg(self):
        """Evaluates system energy and forces using histograms."""
        exec(self.ArgEvalSub.CallObj) 

    def EvalArgStatsThis(self, Weight = 1.):
        """Evaluates statistics of arguments for each energy term;
just for current state."""
        self.Arg.ResetStats()
        exec(self.ArgStatsSub.CallObj)

    def EvalArgStatsAdd(self, Weight = 1.):
        """Evaluates statistics of arguments for each energy term;
adds to running averages."""
        exec(self.ArgStatsSub.CallObj) 

    def EvalArgHistThis(self, Weight = 1.):
        """Makes histograms of arguments for each energy term;
just for current state."""
        self.Arg.ResetHist()
        exec(self.ArgHistSub.CallObj)

    def EvalArgHistAdd(self, Weight = 1.):
        """Makes histograms of arguments for each energy term;
adds to running averages."""
        exec(self.ArgHistSub.CallObj)
        
    def SaveEneState(self, Mode = 0):
        """Saves the current state of force field terms.
Mode = 0 will save for all atoms.
Mode = 1 will save for Sys.TargetAtom
Mode = 2 will save for Sys.TargetMol"""
        exec(self.SaveEneSub.CallObj)
        
    def RevertEneState(self, Mode = 0):
        """Reverts to the old state of force field terms.
Mode = 0 will revert for all atoms.
Mode = 1 will revert for Sys.TargetAtom
Mode = 2 will revert for Sys.TargetMol"""
        exec(self.RevertEneSub.CallObj)

    def EstimateGlobals(self):
        """Estimates force field parameters."""
        Charge = self.Globals.Charge
        NotFixed = np.logical_not(Charge.Fixed)
        n = sum(NotFixed.astype(int))
        if n == 0: return
        Charge[NotFixed] = np.random.rand(n)    
        
    def SetSplineTreatment(self, NonbondCoreEne = None, NonbondCoreEneSlope = None,
                           NonbondCoreMinEne = None, NonbondEneInit = None,
                           BondEneSlope = None, AngleEneSlope = None):
        """Sets defaults in inner and outer regions of spline potentials for
treating high-energy (e.g., core overlap) behavior."""
        EneCore = self.Sys.Units.ConvertEne(NonbondCoreMinEne, self.Sys.TempSet)
        for P in self:
            if not P.IsSpline: continue
            if P.Type == ptypes.PairPotential:
                if P.Filter.Bonded:
                    P.EneSlopeInner = BondEneSlope
                    P.EneSlopeOuter = BondEneSlope
                else:
                    P.EneInner = NonbondCoreEne
                    P.EneOuter = None
                    P.EneSlopeInner = NonbondCoreEneSlope
                    P.EneSlopeOuter = None 
                    P.EneInnerInit = NonbondEneInit
                    P.EneOuterInit = None
                    if NonbondCoreMinEne:
                        P.Knots.Min[0] = EneCore
                        P.Knots.SoftMin[0] = 1.1 * EneCore
                        P.EneInnerInit = 1.1 * EneCore
            elif P.Type == ptypes.AnglePotential:
                P.EneSlopeInner = AngleEneSlope
                P.EneSlopeOuter = AngleEneSlope
    
    def RelaxKnotConstraints(self):
        """Relaxes all inner/outer potential constraints in any spline potential."""
        for P in self:
            if P.IsSpline: P.RelaxKnotConstraints()
             
    def Check(self):
        """Raises an error if any of the cutoffs are greater than half the minimum box length."""
        BoxL = [x for x in self.Sys.BoxL if x > 0.]
        if BoxL:
            HalfL = 0.5 * np.min(BoxL)
            for P in self:
                if not P.Cut is None and P.Cut > HalfL:
                    s = "Cutoff for potential %s is %e which is greater than half the min box length %e." % (P.Label, P.Cut, HalfL)
                    raise ValueError(s)

        
        
