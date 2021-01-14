#/usr/bin/env python

### Class definitions for histograms in SIM suite.
### coded by MSS

import copy
import numpy as np


class Hist1ReportClass(object):
    
    def __init__(self, BinAvgs = None, BinVals = None, Title = "Histogram", 
                 AvgsLabel = "x", ValsLabel = "counts", fmt = r"%11.4e"):
        """Creates a 1D histogram object that can be manipulated or added to other histogram objects.
BinAvgs: list or array of bin averages (e.g., x values)
BinVals: list or array of bin values (e.g., histogram counts or densities)
Title: string title
fmt: default format for histogram entries"""
        self.BinAvgs = []
        self.BinVals = []
        self.AvgsLabel = []
        self.ValsLabel = []
        self.fmt = fmt
        self.Title = Title
        if not BinVals is None and not BinAvgs is None:
            self.AddData(BinVals, ValsLabel, BinAvgs, AvgsLabel)
        
    def AddData(self, BinVals, ValsLabel = "counts", BinAvgs = None, AvgsLabel = None):
        """Adds another column of bin values."""
        if BinAvgs is None:
            if len(self.BinAvgs) == 0:
                raise ValueError("No bin averages have been added to histogram yet.")
            BinAvgs = self.BinAvgs[-1]
        else:
            BinAvgs = copy.deepcopy(BinAvgs)
        if AvgsLabel is None:
            AvgsLabel = self.AvgsLabel[-1]
        if not len(BinAvgs) == len(BinVals):
            raise ValueError("BinAvgs (length %d) and BinVals (length %d) must be same length." % (len(BinAvgs), len(BinVals)))
        BinVals = copy.deepcopy(BinVals)
        self.BinAvgs.append(BinAvgs)
        self.BinVals.append(BinVals)
        self.AvgsLabel.append(AvgsLabel)
        self.ValsLabel.append(ValsLabel)
        
    def __add__(self, other):
        """Adds another histogram class to this one."""
        if not isinstance(other, Hist1ReportClass):
            raise TypeError("Can only add other HistogramReportClasses")
        newself = copy.deepcopy(self)
        for itm in zip(other.BinVals, other.ValsLabel, other.BinAvgs, other.AvgsLabel):
            newself.AddData(*itm)
        return newself
        
    def HistString(self, ShowAllAvgs = False):
        """Returns a string-formatted histogram."""
        l = len(self.fmt % 1)
        s = "%s\n" % self.Title
        NHist = len(self.BinAvgs)
        MaxNBin = max([len(BinAvgs) for BinAvgs in self.BinAvgs])
        #check to see if we need to show bin averages for each histogram
        Show = [True]*NHist
        if not ShowAllAvgs:
            LastBinAvgs = []
            for (i, BinAvgs) in enumerate(self.BinAvgs):
                if not len(BinAvgs) == len(LastBinAvgs):
                    Show[i] = True
                    LastBinAvgs = BinAvgs
                elif np.allclose(LastBinAvgs, BinAvgs):
                    Show[i] = False
                else:
                    Show[i] = True
                    LastBinAvgs = BinAvgs
        if "-" in self.fmt:
            sfmt = "%-" + "%d" % l + "s"
        else:
            sfmt = "%" + "%d" % l + "s"
        for i in range(NHist):
            if Show[i]:
                s += sfmt % self.AvgsLabel[i] + " "
            s += sfmt % self.ValsLabel[i] + " "
        s += "\n"
        for j in range(MaxNBin):
            for i in range(NHist):
                if Show[i]:
                    if j < len(self.BinAvgs[i]):
                        s += self.fmt % self.BinAvgs[i][j] + " "
                    else:
                        s += " "*l + " "
                if j < len(self.BinVals[i]):
                    s += self.fmt % self.BinVals[i][j] + " "
                else:
                    s += " "*l + " "
            s += "\n"
        return s  
        
    def __str__(self):
        return self.HistString()
        
        
        