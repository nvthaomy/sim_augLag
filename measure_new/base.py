#/usr/bin/env python


### Base for measuring objects in SIM suite.
### coded by MSS

import numpy as np

import sys

import basefortran


class MeasureClass(object):

    IsAvg = False
    
    ModuleVars = ["StepFreq", "CycleFreq", "Active", "Val"]
    ModuleConsts = ["NVal" ]
    
    def __init__(self, Sys, Name, StepFreq = 0, CycleFreq = 0, Active = True, 
                 NVal = 1, ValName = None, Report = True, Source = ""):
        """Defines base class for measurements.  Use this class for measurements
that depend on time-dependent properties, like the mean squared displacement.
Use MeasureAvgClass instead for thermodynamic properties / averages.
Name: name of measurement
StepFreq: frequency of updating in terms of integration steps
CycleFreq: frequency of updating in terms of integration cycles
NVal: the number of values given for each measurement
ValName: the name of the value (NVal=1) or a list of names (NVal > 1)
Report: whether or not to include this average in reporting summaries
Active: whether or not measure is active
Source: source code"""
        import copy
        object.__setattr__(self, "LibVars", {})
        self.Sys = Sys
        self.Name = Name
        self.StepFreq = StepFreq
        self.CycleFreq = CycleFreq
        self.Active = Active
        self.NVal = NVal
        self.Val = np.zeros(self.NVal, float)
        self.Avg = np.zeros(self.NVal, float)
        self.ValName = ValName
        self.FortCode = basefortran.NewFortCode(Source)
        self.Report = Report
        self.Filter = None

    def Cleanup(self):
        """Cleans up self."""
        self.LibVars.clear()
        self.Sys = None

    def On(self):
        """Turns measure on."""
        self.Active = True
        self.Init()

    def Off(self):
        """Turns measure off."""
        self.Active = False        
        
    def Init(self):
        """Initializes measure, used before integration."""
        pass
    
    def Finalize(self):
        """Finalizes measure, used after integration."""
        pass    

    def PreLoad(self):
        """Run before compilation."""
        self.Sys.Lib.VarPreLoad(self, self.ModuleVars, self.ModuleConsts, self.FortPrefix)
        
    def PostLoad(self):
        """Run after compilation."""
        self.Sys.Lib.VarPostLoad(self)
        
    def Reset(self):
        """Resets average counters and histogram (if defined)."""
        pass 
    
    def GetNames(self):
        """Returns a list of value names."""
        if self.ValName is None:
            if self.NVal > 1:
                return ["%s_%d" % (self.Name, i) for i in range(self.NVal)]
            else:
                return [self.Name]
        elif type(self.ValName) is list:
            if not len(self.ValName) == self.NVal:
                raise ValueError("ValName is not length NVal")
            return self.ValName
        elif type(self.ValName) is str:
            if self.NVal > 1:
                return ["%s_%d" % (self.ValName, i) for i in range(self.NVal)]
            else:
                return [self.ValName]

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
            

class MeasureAvgClass(MeasureClass):
    
    ModuleVars = MeasureClass.ModuleVars + ["ValSum", "ValSumSq", "Count"]
    ModuleConsts = MeasureClass.ModuleConsts

    IsAvg = True
    
    def __init__(self, Sys, Name, StepFreq = 0, CycleFreq = 0, Active = True,
                 NVal = 1, ValName = None, UseCovar = False,
                 Filter = None, Report = True, Source = ""):
        """Defines base class for measurements that are averages and that do
not depend on time-varying properties.
Name: name of measurement
StepFreq: frequency of updating in terms of integration steps
CycleFreq: frequency of updating in terms of integration cycles
NVal: the number of values given for each measurement, can be matrix
ValName: the name of the value (NVal=1) or a list of names (NVal > 1)
UseCovar: True to compute co-variance between values if NVal > 1
Filter: Filter object specifying to what kinds of atoms to apply this measure
Report: whether or not to include this average in reporting summaries
Active: whether or not measure is active
Source: source code"""
        import copy
        MeasureClass.__init__(self, Sys, Name, StepFreq = StepFreq, CycleFreq = CycleFreq,
                              Active = Active, NVal = NVal, ValName = ValName,
                              Report = Report, Source = Source)      
        del self.Avg
        self.Filter = Filter
        if type(NVal) is tuple:
            if UseCovar:
                raise ValueError("Cannot use covariance if NVal is tuple.")
            self.UseCovar = False
        else:
            self.UseCovar = UseCovar and NVal > 1
        self.Count = 0.
        self.ValSum = np.zeros(self.NVal, float)
        if self.UseCovar:
            self.ValSumSq = np.zeros([self.NVal, self.NVal], float)
        else:
            self.ValSumSq = np.zeros(self.NVal, float)


    def PreLoad(self):
        """Run before compilation."""
        if hasattr(self, "Hist"):
            self.ModuleConsts = self.ModuleConsts + ["HistNBin"]
            self.ModuleVars = self.ModuleVars + ["HistMin", "HistMax", 
                                                 "HistBin", "HistiBin", 
                                                 "HistCount", "Hist"]
        MeasureClass.PreLoad(self)
         

    def SetupHist(self, HistMin, HistMax, HistNBin):
        """Sets up histogram running for this measure;
must be performed before library loading and can only be used
with NVal=1 (1D histogram) or NVal=2 (2D histogram).
HistMin: minimum value for histogram
HistMax: maximum value for histogram
HistNBin: number of histogram bins"""
        if len(self.LibVars):
            raise AttributeError("Cannot modify histogram after library loading.")
        if hasattr(self, "Hist"):
            raise AttributeError("Can only setup histogram once.")
        if not hasattr(self, "NVal"):
            raise AttributeError("Cannot setup histogram for this measure.")
        if self.NVal is tuple:
            raise AttributeError("Cannot use histograms if NVal is tuple.")
        elif self.NVal > 2 or self.NVal == 0:
            raise AttributeError("Can only do 1 and 2-dimensional histograms.")
        self.HistMin = np.zeros(self.NVal, float)
        self.HistMin[:] = HistMin
        self.HistMax = np.zeros(self.NVal, float)
        self.HistMax[:] = HistMax
        self.HistNBin = np.zeros(self.NVal, int)
        self.HistNBin[:] = HistNBin
        self.Hist = np.zeros(self.HistNBin, float)
        self.HistBin = np.zeros(self.NVal, float)
        self.HistiBin = np.zeros(self.NVal, float)
        self.HistCount = 0.
        self.HistBin[:] = (self.HistMax - self.HistMin) / float(self.HistNBin + 1.e-300)
        self.HistiBin[:] = 1./self.HistBin

    def HistString(self, fmt = "%-11.4e"):
        """Returns a string with the contents of the histogram."""
        if not hasattr(self, "Hist"):
            raise AttributeError("Histogram not defined for this measure.")
        n = len(fmt % 1)    
        Names = self.GetNames()
        if self.NVal == 1:
            BinAvg0 = self.HistMin[0] + self.HistBin[0] * np.arange(0.5, self.HistNBin[0], 1.)
            s = Names[0][:n].ljust(n) + "  " + "Freq".ljust(n) + "\n"
            for (a, v) in zip(BinAvg0, self.Hist):
                s += (fmt % a + "  " + fmt % v + "\n")
        else:
            BinAvg0 = self.HistMin[0] + self.HistBin[0] * np.arange(0.5, self.HistNBin[0], 1.)
            BinAvg1 = self.HistMin[1] + self.HistBin[1] * np.arange(0.5, self.HistNBin[1], 1.)
            s = Names[0][:n].ljust(n) + "  " + Names[1][:n].ljust(n) + "  ""Freq".ljust(n) + "\n"
            for i in range(self.HistNBin[0]):
                for j in range(self.HistNBin[1]):
                    s += (fmt % BinAvg0[i] + "  " + fmt % BinAvg1[j]
                          + "  " + fmt % self.Hist[i,j] + "\n")
        return s

    def ResetHist(self):
        """Resets histogram running."""
        if hasattr(self, "Hist"):
            self.Hist.fill(0.)
            self.HistCount = 0.

    def ResetAvgs(self):
        """Resets average counters to zero."""
        self.ValSum.fill(0.)
        self.ValSumSq.fill(0.)
        self.Count = 0.

    def Reset(self):
        """Resets average counters and histogram (if defined)."""
        self.ResetAvgs()
        if hasattr(self, "Hist"):
            self.ResetHist()      

    def __getattr__(self, name):
        Var = self.LibVars.get(name, None)
        if not Var is None:
            if not Var.shape: Var = Var.item()
            return Var
        elif name == "Avg":
            ret = self.ValSum / (self.Count + 1.e-300)
            if self.NVal == 1:
                return ret[0]
            else:
                return ret
        elif name == "AvgSq":
            ret = self.ValSumSq / (self.Count + 1.e-300)
            if self.NVal == 1:
                return ret[0]
            else:
                return ret
        elif name == "Var":
            if self.UseCovar:
                ret = self.ValSumSq.diagonal() / (self.Count + 1.e-300) \
                      - (self.ValSum / (self.Count + 1.e-300))**2
            else:
                ret = self.ValSumSq / (self.Count + 1.e-300) \
                      - (self.ValSum / (self.Count + 1.e-300))**2
            if self.NVal == 1:
                return ret[0]
            else:
                return ret
        elif name == "Std":
            if self.UseCovar:
                ret = self.ValSumSq.diagonal() / (self.Count + 1.e-300) \
                      - (self.ValSum / (self.Count + 1.e-300))**2
            else:
                ret = self.ValSumSq / (self.Count + 1.e-300) \
                      - (self.ValSum / (self.Count + 1.e-300))**2
            ret = np.sqrt(ret)
            if self.NVal == 1:
                return ret[0]
            else:
                return ret
        elif name == "Covar":
            if not self.UseCovar:
                raise AttributeError("Covar not defined for this measure.")
            ret = self.ValSumSq / (self.Count + 1.e-300) \
                  - np.outer(self.ValSum, self.ValSum) / (self.Count + 1.e-300)**2
            return ret
        else:
            raise AttributeError(name)

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
            if name in ["HistMin", "HistMax"]:
                self.HistBin[:] = (self.HistMax - self.HistMin) / float(self.HistNBin + 1.e-300)
                self.HistiBin[:] = 1./self.HistBin



def ListMeasures(obj):
    """Returns a list of all measures in obj (e.g., a module)"""
    l = []
    if type(obj) is dict:
        Vars = obj.values()
    else:
        Vars = [getattr(obj, Name, None) for Name in dir(obj)]
    for Var in Vars:
        try:
            if issubclass(Var, MeasureClass):
                l.append(Var)
        except TypeError:
            pass
    return l


class MeasuresClass(list):
    
    def __init__(self, Sys):
        """Holder object for measurements.
Sys: System object to which to apply measurements"""
        list.__init__(self)
        self.Sys = Sys
        #fortran subroutine for energy
        from sim.fortran import Subroutine
        self.Sub = Subroutine("calcmeasures", SubStr = "self.Sub", PreCallCode = "Sys = self.Sys")
        self.Sub.PreLoad = self.PreLoad
        self.Sub.PostLoad = self.PostLoad
        self.Sys.Lib.append(self.Sub)
        self.fmt = "%11.4e"
        self.IntAction = None
        
    def __getattr__(self, name):
        for M in self:
            if M.Name == name:
                return M
        raise AttributeError(name)
        
    def Show(self):
        """Shows a list of measures."""
        print "Available measurements:"
        print "  Currently ACTIVE"
        for M in self:
            if M.Active:
                print "    %s" % M.Name
        print "  Currently INACTIVE"
        for M in self:
            if not M.Active:
                print "    %s" % M.Name

    def Cleanup(self):
        for M in self:
            M.Cleanup()
        self.Sys = None
        self.Sub = None

    def AllOn(self):
        """Turns all measures on."""
        for M in self:
            M.On()

    def AllOff(self):
        """Turns all measures off."""
        for M in self:
            M.Off()
            
    def Init(self):
        """Initializes measures, used before integration."""
        for M in self:
            if M.Active: M.Init()

    def Finalize(self):
        """Finalizes measures, used after integration."""
        for M in self:
            if M.Active: M.Finalize() 
    
    def Reset(self):
        """Resets all counters in all measures."""
        for M in self:
            M.Reset()

    def ResetAvgs(self):
        """Resets averages in all measures."""
        for M in self:
            M.ResetAvgs()

    def ResetHist(self):
        """Resets histogram in all measures."""
        for M in self:
            M.ResetHist()            

    def PreLoad(self):
        """Prepares class for compiling."""
        #run pre-load command for all measures
        for (ind, M) in enumerate(self):
            M.FortPrefix = "M%d_" % ind
            M.PreLoad()      
        #add the number of measures
        self.Sys.Lib.AddVars("int NMeasure = %d" % len(self), IsConstant = True)
        #main source for energy routine
        fc = basefortran.GetFortCode(self, self.Sys)
        self.Sub.AddFortCode(fc)   

    def PostLoad(self):
        """Run after compilation."""
        for M in self:
            M.PostLoad()

    def Eval(self, StepNum = 0, CycleNum = 0, MeasureAll = False, Weight = 1.):
        """Runs all measurements.  Called from an integrator.
StepNum: current step number
CycleNum: current cycle number
MeasureAll: override individual UseMeasure for each measurement
            and measure all of them for this call
Weight: weight of this evaluation in averages"""
        exec(self.Sub.CallObj)

    def ReportHead(self, MaxWidth = None):
        """Returns a header string for reporting of averages.
MaxWidth: maximum width of string in number of characters"""
        l = len(self.fmt % -1)
        s = []
        for M in self:
            if not M.Active or not M.Report: continue
            Names = M.GetNames()
            Names = [nm.rjust(l)[:l] for nm in Names]
            s.extend(Names)
        s = " ".join(s)
        if not MaxWidth is None: s = s[:MaxWidth]
        return s

    def ReportAvgs(self, MaxWidth = None):
        """Returns a string of reporteds averages.
MaxWidth: maximum width of string in number of characters"""
        s = []
        for M in self:
            if not M.Active or not M.Report: continue
            if M.NVal == 1:
                s.append(self.fmt % M.Avg)
            else:
                Avg = M.Avg
                for i in range(M.NVal):
                    s.append(self.fmt % Avg[i])
        s = " ".join(s)
        if not MaxWidth is None: s = s[:MaxWidth]
        return s
    
    def ReportVals(self, MaxWidth = None):
        """Returns a string of instantaneous values.
MaxWidth: maximum width of string in number of characters"""
        s = []
        for M in self:
            if not M.Active or not M.Report: continue
            for i in range(M.NVal):
                s.append(self.fmt % M.Val[i])
        s = " ".join(s)
        if not MaxWidth is None: s = s[:MaxWidth]
        return s 
        
    def VerboseOutput(self, StepFreq = 0, CycleFreq = 0, Int = None, fobj = sys.stdout):
        """Adds periodic output to current integrator."""
        from sim.integrate import Action
        if Int is None:
            Int = self.Sys.Int
        if not self.IntAction is None:
            Int.DelAction(self.IntAction)
            self.IntAction = None
        if StepFreq > 0 or CycleFreq > 0:
            #add an action to happen during integration
            def PrintHead(Sys):
                Sys.Measures.Reset()
                fobj.write("%7s " % "Step" + Sys.Measures.ReportHead() + "\n")
            def PrintMeasures(Sys):
                fobj.write("%7d " % Sys.Int.StepNum + Sys.Measures.ReportVals() + "\n")
            def PrintAverages(Sys):
                s = "%7s " % "avgs" + Sys.Measures.ReportAvgs()
                n = len(s)
                s = "-"*n + "\n" + s + "\n" + "="*n + "\n"
                fobj.write(s)
            self.IntAction = Action(StepFreq = StepFreq, CycleFreq = CycleFreq,
                                    Fn = PrintMeasures, Name = "MeasuresOutput",
                                    InitFn = PrintHead, FinalFn = PrintAverages)
            Int.Actions.append(self.IntAction)