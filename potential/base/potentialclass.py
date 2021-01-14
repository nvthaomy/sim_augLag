### Class definitions for potential energy functions in SIM suite.
### coded by MSS

import numpy as np

import potentialarray as PA
import argarray as AA

import sim.fortran as fortran
import potentialtypes as ptypes


class PotentialClass(object):
    """Generic class for developing potential energy functions."""
    
    Names = ["base"]
    Type = None
    #is this a spline?  by default, off
    IsSpline = False
    #does this potential use atom charges?
    UsesCharge = False
    #does this potential require atom types outside of filter (e.g., charges)?
    UsesATypes = False
    #can argument histograms not be used to evaluate potential?
    NoHistEval = False
    #don't automatically scale energies/virial in updating potential for molecule delete/ins?
    NoEnergyUpdateScaling = False
    #whether or not to use energy fluctuations (Experimental)
    UsesFluct = False    


    def __init__(self, Sys, Label = "", Cut = None,
                 Filter = None, NType = 1, TypeLabels = None, Source = None):
        object.__setattr__(self, 'LibVars', {})
        self.Sys = Sys
        #parameters for this potential
        self.Param = PA.Parray(Name = Label)
        #argument stats for this potential
        self.Arg = AA.ArgArray(NType = NType, TypeLabels = TypeLabels)
        #label
        self.Label = Label
        #fortran code
        #Fortran code object
        self.FortCode = fortran.FortCode()        
        if not Source is None:
            self.SetSource(Source)
        #cutoff
        self.Cut = Cut
        #set the filter
        self.Filter = Filter
        #set default reporting
        self.SetReport()
        #set the list of frozen parameters
        self.__FrozenInd = []
        #set the estimated histogram fractional error
        self.HistFracErr = 1.0
        #fortran module vars
        self.ModuleVars = []
        self.ModuleConsts = []
        
    def SetSource(self, Source):
        """Adds source code to this potential."""
        self.FortCode.Add(Source = Source)

    def SetReport(self, ArgMin = None, ArgMax = None, ArgDelta = None, ArgNBin = 100):
        """Updates settings for reporting potential summaries."""
        if ArgMin is None:
            if self.Type == ptypes.TorsionPotential:
                self.Arg.ReportMin = -np.pi
            else:
                self.Arg.ReportMin = 0.
        else:
            self.Arg.ReportMin = ArgMin
        if ArgMax is None:
            if self.Type == ptypes.PairPotential and not self.Cut is None:
                self.Arg.ReportMax = self.Cut
            elif self.Type in [None, ptypes.FieldPotential, ptypes.PairPotential]:
                self.Arg.ReportMax = 10. * self.Sys.Units.LScale
            elif self.Type in [ptypes.AnglePotential, ptypes.TorsionPotential]:
                self.Arg.ReportMax = np.pi
            elif self.Type == ptypes.LocalDensityPotential:
                self.Arg.ReportMax = 20.
            elif self.Type == ptypes.GlobalPotential:
                self.Arg.ReportMax = 10.
        else:
            self.Arg.ReportMax = ArgMax
        if ArgDelta is None:
            self.Arg.ReportDelta = (self.Arg.ReportMax - self.Arg.ReportMin) / ArgNBin
        else:
            self.Arg.ReportDelta = ArgDelta

    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        self.Param = None
        self.Sys = None        

    def __setattr__(self, name, val):
        if name == "Filter" and not val is None:
            Filter = val
            if self.Type == ptypes.FieldPotential:
                if not Filter.N == 1:
                    raise TypeError("Filter must be for single atoms.")
            elif self.Type == ptypes.PairPotential:
                if not Filter.N == 2:
                    raise TypeError("Filter must be for atom pairs.")
                elif not Filter.MaxBondOrd <= 4:
                    raise TypeError("Cannot parse bond orders beyond 4.")
            elif self.Type == ptypes.AnglePotential:
                if not (Filter.N == 3 and Filter.Bonded == True):
                    raise TypeError("Filter must be for bonded atom triples.")
                elif not Filter.MaxBondOrd is None:
                    raise TypeError("Cannot parse bond orders for angle potentials.")
            elif self.Type == ptypes.TorsionPotential:
                if not (Filter.N == 4 and Filter.Bonded == True):
                    raise TypeError("Filter must be for bonded atom quartets.")
                elif not Filter.MaxBondOrd is None:
                    raise TypeError("Cannot parse bond orders for torsion potentials.")
            elif self.Type == ptypes.LocalDensityPotential: 
                if not Filter.N == 2:
                    raise TypeError("Filter must be for atom pairs.")
                if not Filter.Ordered:
                    raise TypeError("LocalDensityPotential must have Filter.Ordered = True")
            elif self.Type == ptypes.GlobalPotential:
                pass
            else:
                raise TypeError("I don't recognize this kind of potential.")
            object.__setattr__(self, "Filter", Filter)
        elif name in self.__dict__:
            object.__setattr__(self, name, val)        
        else:
            #check for a module var or parameter
            Var = self.LibVars.get(name, None)
            if "Param" in self.__dict__:
                IsParam = hasattr(self.Param, name)
            else:   
                IsParam = False
            if IsParam:
                #it is a parameter
                setattr(self.Param, name, val)
            elif Var is None:
                #it is a normal variable
                v = self.__dict__.get(name, None)
                if isinstance(v, np.ndarray):
                    v[:] = val
                else:
                    object.__setattr__(self, name, val) 
            else:
                #it is a module variable
                if Var.shape:
                    Var[:] = val
                else:
                    Var.itemset(val)

    def __getattr__(self, name):
        if hasattr(self.Param, name):
            return getattr(self.Param, name)
        else:
            Var = self.LibVars.get(name, None)
            if Var is None:
                raise AttributeError(name) 
            else:
                if not Var.shape: Var = Var.item()
                return Var

    def __repr__(self):
        if self.Label:
            return "<%s>" % self.Label
        else:
            return repr(self.__class__)

    def PreLoad(self):
        """Actions to take before loading/compiling Fortran lib."""
        self.Sys.Lib.VarPreLoad(self, self.ModuleVars, self.ModuleConsts, self.FortPrefix)

    def PostLoad(self):
        """Actions to take after loading/compiling Fortran lib."""
        self.Sys.Lib.VarPostLoad(self)

    def SetParam(self):
        """Function for setting parameters."""
        pass

    def Update(self):
        """Updates potential given parameter changes."""
        pass

    def SetupTest(self):
        """Sets parameters for a test run."""
        pass

    def Estimate(self, Beta):
        """Estimates parameters from arg count."""
        pass

    def SetBounds(self, MaxFracChange = 0.1):
        """Sets bounds on parameters and parameter changes based on argument ranges and scales."""
        pass
    
    def GetTypeInd(self, *args, **kwargs):
        """Returns the type number for this potential."""
        #to be defined by subclasses
        return 0 
        
    def SetTypeInd(self, TypeInd):
        """Sets the type to be used in Val and DVal, or None for generic."""
        pass
    
    def Val(self, x):
        """Returns the value of the potential."""
        #to be defined by subclasses
        return 0
        
    def DVal(self, x):
        """Returns the derivative of the potential."""
        #to be defined by subclasses
        return 0
    
    def FLVal(self,x):
        """Returns the fluctuation value of the potential."""
        if not self.UsesFluct:
            raise ValueError("Potential does not have a fluctuation term.")
        #remainder to be defined by subclass
        return 0

    def ReportString(self, fmt = r"%12.4e"):
        """Returns a string of arguments and values."""
        l = len(fmt % 1)
        s = ""
        s += "POTENTIAL %s\n" % self.Name
        if self.UsesFluct:
            s += "%-*s %-*s %-*s %-*s\n" % (l, 'arg', l, 'val', l, 'deriv', l, 'fluct')
        else:
            s += "%-*s %-*s %-*s\n" % (l, 'arg', l, 'val', l, 'deriv')
        ReportMin = self.Arg.ReportMin.min()
        ReportMax = self.Arg.ReportMax.max()
        ReportDelta = self.Arg.ReportDelta.min()
        for x in np.arange(ReportMin, ReportMax, ReportDelta):
            if self.UsesFluct:
                s += fmt % x + " " + fmt % self.Val(x) + " " + fmt % self.DVal(x) + " " + fmt % self.FLVal(x) +  '\n'
            else:
                s += fmt % x + " " + fmt % self.Val(x) + " " + fmt % self.DVal(x) + '\n'
        return s

    def HistString(self, fmt = r"%-11.4e", ArgLabel = 'arg', HistLabel = 'hist'):
        """Returns a string of arguments histograms."""
        return self.Arg.HistString(Name = "POTENTIAL %s" % self.Name, fmt=fmt,
                                   ArgLabel = ArgLabel, HistLabel = HistLabel)
        
    def HistList(self, fmt = r"%-11.4e", ArgLabel = 'arg', HistLabel = 'hist'):
        """Returns a list of arguments histograms."""
        return self.Arg.HistList(Name = "POTENTIAL %s" % self.Name, fmt=fmt,
                                 ArgLabel = ArgLabel, HistLabel = HistLabel)
      
    def ParamDict(self):
        """Returns a dictionary of parameters."""
        d = {}
        for key in self.Param.Maps:
            val = getattr(self, key)
            if len(val) > 1:
                val = list(val)
            else:
                val = float(val)
            d[key] = val
        return d                
      
    def ParamString(self, fmt = r"%-11.4e"):
        """Returns a string summarizing parameters."""
        d = self.ParamDict()
        l = []
        for (k,v) in d.items():
            if type(v) is float:
                s = "'%s' : " % k + fmt % v
            else:
                s = "'%s' : [" % k + ", ".join([fmt % x for x in v]) + "]"
            if len(l):
                l.append(" " + s)
            else:
                l.append("{" + s) 
        if len(l):
            return ",\n".join(l) + "}"
        else:
            return "{}"
    def FreezeSpecificParam(self,paramFix):
        """Freezes specified parameters temporarily."""
        # - Modification by Nick Sherck
        # - i runs from 0 to # parameters in potential object
        # - User specifies paramFix which is a list of indexes specifiy which of the 
        #   parameters in i to fix. (not the names)
        for (i,Fixed) in enumerate(self.Param.Fixed):
            for j in paramFix:
                if i==j:
                    print "i j"
                    print i,j
                    self.__FrozenInd.append(i)
                    self.Param.Fixed[i] = True
                else:
                    pass


    def FreezeParam(self):
        """Freezes any unfixed parameters temporarily."""
        for (i,Fixed) in enumerate(self.Param.Fixed):
            if not Fixed:
                self.__FrozenInd.append(i)
                self.Param.Fixed[i] = True

    def FreezeSpecificParam(self,paramFix):
        """Freezes specified parameters temporarily."""
        # - i runs from 0 to # parameters in potential object
        # - User specifies paramFix which is a list of indexes specifying which of the 
        #   parameters (in self.Param) to fix. (not the names)
        #   the order of parameters in self.Param is the order in which they were added
        #   in the force field __init__ function
        print "i,j"
        for (i,Fixed) in enumerate(self.Param.Fixed):
            for j in paramFix:
                if i==j:
                    print i, j
                    self.__FrozenInd.append(i)
                    self.Param.Fixed[i] = True
                else:
                    pass
    
    def UnfreezeParam(self):
        """Returns to normal any temporarily frozen parameters."""
        while len(self.__FrozenInd):
            self.Param.Fixed[self.__FrozenInd.pop(0)] = False     
