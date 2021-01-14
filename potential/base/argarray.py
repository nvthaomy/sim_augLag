#/usr/bin/env python

### Class definitions for argument arrays in SIM suite.
### coded by MSS

import numpy as np
from sim.histogram import Hist1ReportClass

SetupHistBuffer = 1.e-6


def npStrides(Strides, dtype):
    itemsize = np.dtype(dtype).itemsize
    return tuple([itemsize*x for x in Strides])

def Strides(a):
    return tuple([x/a.itemsize for x in a.strides])

def LinkArray(Array, Master, StartInd):
    n = Array.size
    if Array.flags.c_contiguous:
        Master[StartInd:StartInd + n] = Array.flatten('C')
    else:
        Master[StartInd:StartInd + n] = Array.flatten('F')
    Array.data = Master[StartInd:StartInd + n].data

def MakeMasterArray(ArrayList):
    n = sum([a.size for a in ArrayList])
    dtype = ArrayList[0].dtype
    Master = np.zeros(n, dtype=dtype)
    Ind = 0
    for a in ArrayList:
        LinkArray(a, Master, Ind)
        Ind += a.size
    return Master

   

class ArgArray(object):
    """Generic class for describing and containing argument stats."""
    Array1Vars = ["Min", "Max", "Count", "Sum", "SumSq",
                  "ReportMin", "ReportMax", "ReportDelta",
                  "HistMin", "HistMax", "HistBinw", "HistiBinw"]
    Array2Vars = ["Hist"]
    SettingsVars = ["HistNBin", "HistMin", "HistMax", "HistReportNBin"]
    ModuleVars = ["WeightSumHist", "WeightSumStats"] + Array1Vars + Array2Vars
    ModuleConsts = ["HistNBin"]

    def __init__(self, NType = 1, HistNBin = 0, TypeLabels = None):
        object.__setattr__(self, "LibVars", {})
        self.Ind = None
        self.Ind2 = None
        self.HistNBin = HistNBin
        self.HistReportNBin = None
        self.NType = NType
        if TypeLabels is None:
            self.TypeLabels = ["Type%d" % i for i in range(NType)]
        elif len(TypeLabels) != NType:
            raise ValueError("TypeLabels does not have NType elements")
        else:
            self.TypeLabels = TypeLabels
        for v in ArgArray.Array1Vars:
            setattr(self, v, np.zeros(self.NType, dtype=float))
        for v in ArgArray.Array2Vars:
            setattr(self, v, np.zeros((self.NType,0), dtype=float, order = 'C'))
        self.N = self.Min.size
        self.N2 = self.Hist.size
        self.Reset()
        
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
            if name == "Avg":
                return self.Sum / (self.Count + 1.e-6)
            elif name == "Std":
                Avg = self.Sum / (self.Count + 1.e-6) 
                AvgSq = self.SumSq / (self.Count + 1.e-6)
                return np.sqrt(np.abs(AvgSq - Avg*Avg))
            else:
                raise AttributeError(name) 
        else:
            if not Var.shape: Var = Var.item()
            return Var            

    def PreLoad(self, Lib):
        """Runs before compilation."""
        Lib.VarPreLoad(self, self.ModuleVars, self.ModuleConsts, self.FortPrefix)
        
    def PostLoad(self, Lib):
        Lib.VarPostLoad(self)        

    def SetupHist(self, NBin = None, Min = None, Max = None, ReportNBin = None):
        """Sets up a histogram for this argument."""
        if not NBin is None:
            if NBin != self.HistNBin and not self.Ind is None:
                raise AttributeError("Cannot modify histogram size after library loading.")
            else:
                del self.HistNBin
                self.HistNBin = NBin
                del self.Hist
                self.Hist = np.zeros((self.NType,self.HistNBin),
                                     dtype=float, order = 'C')
                self.N2 = self.Hist.size
        if not ReportNBin is None:
            if ReportNBin <= 0:
                self.HistReportNBin = None
            else:
                if self.HistNBin % ReportNBin > 0:
                    raise ValueError("ReportNBin must be evenly divisible into NBin")
                self.HistReportNBin = ReportNBin
        if Min is None:
            self.HistMin[:] = self.Min[:] - SetupHistBuffer * np.abs(self.Min[:])
        else:
            self.HistMin[:] = Min - SetupHistBuffer * np.abs(Min)
        if Max is None:
            self.HistMax[:] = self.Max[:] + SetupHistBuffer * np.abs(self.Max[:])
        else:
            self.HistMax[:] = Max + SetupHistBuffer * np.abs(Max)
        if self.HistNBin == 0:
            self.HistBinw[:] = 0.
            self.HistiBinw[:] = 0.
        else:
            self.HistBinw[:] = (self.HistMax - self.HistMin) / self.HistNBin
            self.HistiBinw[:] = np.where(self.HistBinw == 0, 0, 1. / self.HistBinw)

    def ResetHist(self):
        """Resets the histogram."""
        self.Hist.fill(0.)
        self.WeightSumHist = 0.

    def ResetStats(self):
        """Resets statistics."""
        self.Min.fill(1.e300)
        self.Max.fill(-1.e300)
        self.Sum.fill(0.)
        self.SumSq.fill(0.)
        self.Count.fill(0.)
        self.WeightSumStats = 0.

    def Reset(self):
        """Resets histogram and statistics."""
        self.ResetHist()
        self.ResetStats()
        
    def NormalizeHist(self):
        """Normalizes the histogram based on added entries and weights."""
        self.Hist /= self.WeightSumHist
        self.WeightSumHist = 1.

    def NormalizeStats(self):
        """Normalizes statistics based on added entries and weights."""
        self.Sum /= self.WeightSumStats
        self.SumSq /= self.WeightSumStats
        self.Count /= self.WeightSumStats
        self.WeightSumStats = 1.
        for (i, (Min, Max)) in enumerate(zip(self.Min, self.Max)):
            if Min > Max:
                self.Min[i] = 1.
                self.Max[i] = 1.

    def Normalize(self):
        """Normalizes histogram and statistics based on added entries and weights."""
        self.NormalizeStats()
        self.NormalizeHist()
        
        
    def HistList(self, Name = 'ARGUMENTS', fmt = r"%-11.4e", 
                 ArgLabel = 'arg', HistLabel = 'hist', NormalizeDensity = True):
        """Returns a list of histogram report objects."""
        l = []
        for i in range(self.NType):
            if self.NType == 1:
                Title = "%s\n" % Name
            else:
                Title = "%s(%s)\n" % (Name, self.TypeLabels[i])
            if self.HistReportNBin is None:
                BinAvgs = self.HistMin[i] + self.HistBinw[i] * np.arange(0.5, self.HistNBin, 1.0)
                BinVals = self.Hist[i,:]
            else:
                if self.HistNBin % self.HistReportNBin > 0:
                    raise ValueError("HistReportNBin must be evenly divisible into NBin")
                r = self.HistNBin / self.HistReportNBin
                BinAvgs = self.HistMin[i] + self.HistBinw[i] * r * np.arange(0.5, self.HistReportNBin, 1.0)
                BinVals = self.Hist[i,:].copy().reshape((-1,r), order='C').sum(axis=1)
            if NormalizeDensity and len(BinVals) > 1:
                BinVals = BinVals.copy() / (BinAvgs[1] - BinAvgs[0])
            h = Hist1ReportClass(BinAvgs = BinAvgs, BinVals = BinVals, 
                                 Title = Title, fmt = fmt,
                                 AvgsLabel = ArgLabel, ValsLabel = HistLabel)
            l.append(h)
        return l          
    
    def HistString(self, *arg, **kwarg):
        """Returns a string of argument histograms."""
        return "\n".join([h.HistString() for h in self.HistList(*arg, **kwarg)]) 
           

            
class MasterArgArray(object):
    """Generic master class for describing and containing argument stats."""
    
    def __init__(self, ArgArrayList):
        self.ArgArrayList = ArgArrayList

    def NormalizeHist(self):
        """Normalizes the histogram based on added entries and weights."""
        for a in self.ArgArrayList:
            a.NormalizeHist()

    def ResetHist(self):
        """Resets the histogram."""
        for a in self.ArgArrayList:
            a.ResetHist()

    def NormalizeStats(self):
        """Normalizes statistics based on added entries and weights."""
        for a in self.ArgArrayList:
            a.NormalizeStats()

    def ResetStats(self):
        """Resets statistics."""
        for a in self.ArgArrayList:
            a.ResetStats()

    def Normalize(self):
        """Normalizes histogram and statistics based on added entries and weights."""
        for a in self.ArgArrayList:
            a.Normalize()

    def Reset(self):
        """Resets histogram and statistics."""
        for a in self.ArgArrayList:
            a.Reset()
    
    def GetArgRanges(self):
        """Gets the histogram mins, maxs, counts for all arguments."""
        Mins = [a.Min.copy() for a in self.ArgArrayList]
        Maxs = [a.Max.copy() for a in self.ArgArrayList]
        return Mins, Maxs      
    
    def SetHistRanges(self, Mins, Maxs):
        """Sets the histogram mins, maxs for all arguments."""
        for (a, Min, Max) in zip(self.ArgArrayList, Mins, Maxs):
            a.SetupHist(Min = Min , Max = Max)      
            
    def GetData(self):
        """Gives current histogram data for later use."""
        import copy
        ret = []
        for a in self.ArgArrayList:
            d = dict([[var, copy.copy(getattr(a, var))] for var in ArgArray.ModuleVars])
            ret.append(d)
        return ret
        
    def SetData(self, HistData):
        """Sets current histogram data using saved tuple."""
        for (i, a) in enumerate(self.ArgArrayList):
            for (k, v) in HistData[i].iteritems():
                setattr(a, k, v)


def AddToData(HistDataTot, HistData, Weight = 1.):
    """Adds all data to a master data object, with weight as indicated.
This is used for averaging histograms of different systems.
Returns updated HistDataTot."""
    if HistDataTot is None:
        HistDataTot = []
        for itm in HistData:
            d = {}    
            for (k,v) in itm.iteritems():
                d[k] = v * Weight
            HistDataTot.append(d)
    else:
        for (d, itm) in zip(HistDataTot, HistData):
            for (k,v) in itm.iteritems():
                d[k] = d[k] + v * Weight
    return HistDataTot

    