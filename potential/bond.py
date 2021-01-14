#/usr/bin/env python


import base
import numpy as np


class Bond(base.PotentialClass):
    """Harmonic bond interactions."""

    Names = ["bond"]    

    Source = base.basefortran.SourcePairTemplate + """
>>> defs
float val1

>>> pairenergy
val1 = DIJ - Dist0(0)
THISU = FConst(0) * val1*val1

>>> pairvirial
THISW = 2.d0 * FConst(0) * val1*DIJ

>>> pairduparam
DU_FConst(0) = val1 * val1
DU_Dist0(0) = -2.d0 * FConst(0) * val1
DDU_Dist0_Dist0(0,0) = 2.d0 * FConst(0)
DDU_Dist0_FConst(0,0) = -2.d0 * val1

>>> pairdwparam
DW_FConst(0) = 2.d0 * val1 * DIJ
DW_Dist0(0) = -2.d0 * FConst(0) * DIJ
DDW_Dist0_FConst(0,0) = -2.d0 * DIJ
"""

    TestArgs = {"Dist0":1.2, "FConst":2.3}
    Type = base.ptypes.PairPotential

    def __init__(self, Sys, Label = "", Filter = None,
                 Fixed = False,
                 Dist0 = None, FConst = None):
        """Initializes a harmonic bond interaction."""
        base.PotentialClass.__init__(self, Sys, Label, Filter = Filter,
                                     Source = Bond.Source)
        if Dist0 is None:
            Dist0 = Sys.Units.LScale
        if FConst is None:
            FConst = Sys.Units.EScale / Sys.Units.LScale**2
        self.Param.Add("Dist0", 1, Value = Dist0, Fixed = Fixed,
                       Scale = self.Sys.Units.LScale)
        self.Param.Add("FConst", 1, Value = FConst, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.EScale / self.Sys.Units.LScale**2)

        self.Update()

    def SetParam(self, Dist0, FConst):
        """Sets parameters for this potential."""
        self.Dist0 = Dist0
        self.FConst = FConst
        self.Update()
        
    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.FConst.SoftMin = 0.0001 * self.Sys.Units.EScale / self.Sys.Units.LScale**2
        if not MaxFracChange is None:
            self.FConst.MaxChange = MaxFracChange * self.Sys.Units.EScale / self.Sys.Units.LScale**2
            self.Dist0.MaxChange = MaxFracChange * self.Sys.Units.LScale

    def Val(self, x):
        """Returns the value of the potential."""
        d = x - self.Dist0[0]
        return self.FConst[0] * d**2
        
    def DVal(self, x):
        """Returns the derivative of the potential."""
        d = x - self.Dist0[0]
        return 2. * self.FConst[0] * d        

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        minx = self.Arg.HistMin
        maxx = self.Arg.HistMax
        dx = self.Arg.HistBinw
        xvals = np.arange(minx + 0.5*dx, maxx, dx)
        hist = self.Arg.Hist[0] / np.sum(self.Arg.Hist[0])
        avgx = sum([h*x for (h,x) in zip(hist, xvals)])
        stdx = sum([h * (x-avgx)**2 for (h,x) in zip(hist, xvals)])  
        k = self.Sys.Units.kB * self.Sys.TempSet / (2. * stdx)
        self.SetParam(avgx, k)


    def DepricEstimate(self):
        """Estimates parameters based on argument ranges and scales."""
        print self.Arg
        minhistval=min([(-np.log(self.Arg.Hist[0][i]/self.Arg.Hist.sum())) for i,x in enumerate(np.linspace(self.Arg.HistMin,self.Arg.HistMax,self.Arg.HistNBin)) if self.Arg.Hist[0][i]>0])
        rval=[round((-np.log(self.Arg.Hist[0][i]/self.Arg.Hist.sum())-minhistval)/((x-self.Arg.Avg)**2)[0]) for i,x in enumerate(np.linspace(self.Arg.HistMin,self.Arg.HistMax,self.Arg.HistNBin)) if self.Arg.Hist[0][i]>0]
        import scipy.stats
        FConst=scipy.stats.mode(rval)[0][0]
        Dist0 = self.Arg.Avg
        print FConst#FConst = self.Sys.Units.EScale / self.Sys.Units.LScale**2
        self.SetParam(Dist0, FConst)

    #def Estimate(self):
    #    """Estimates parameters based on argument ranges and scales."""
    #    Dist0 = self.Arg.Avg
    #    FConst = self.Sys.Units.EScale / self.Sys.Units.LScale**2
    #    self.SetParam(Dist0, FConst)
   


                
            
        
        

        
