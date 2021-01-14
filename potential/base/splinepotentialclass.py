### Class definitions for spline potential energy functions in SIM suite.
### coded by MSS


import numpy as np
from potentialclass import PotentialClass


class SplinePotentialClass(PotentialClass):
    """Generic class for developing spline-based potential energy functions."""
    
    #bspline matrix
    BFunc = np.array([[ 1.,  4.,  1.,  0.],
                      [-3.,  0.,  3.,  0.],
                      [ 3., -6.,  3.,  0.],
                      [-1.,  3., -3.,  1.]], float) / 6.
    
    
    SourceAcyclic = """
>>> defs
float SPt
int SPInd
float SPx
float SPdm1
float SPd0
float SPdp1
float SPdp2
float Thisx
float SPd(4)
float SPThisU
float SPThisW
float SPWCoef

>>> spindex
SPx = (Thisx - SPMin) * SPiDist
SPind = max(min(int(SPx), NInterval - 1), 0) 
SPt = max(min(SPx - real(SPind), 1.d0), 0.d0)

>>> spenergy
SPThisU = SPC0(SPInd) + SPt * (SPC1(SPInd) + SPt * (SPC2(SPInd) + SPt * SPC3(SPInd)))
THISU = SPCoef * SPThisU

>>> spvirial
if (Thisx < SPMin .or. Thisx > SPMax) then
    SPThisW = 0.
else
    SPThisW = (SPC1(SPInd) + SPt * (2.*SPC2(SPInd) + SPt * 3.*SPC3(SPInd))) * SPx
endif
THISW = SPCoef * SPThisW

>>> spdenergy
if (Thisx < SPMin .or. Thisx > SPMax) then
    ThisDU = 0.
else
    ThisdU = SPCoef * (SPC1(SPInd) + SPt * (2.*SPC2(SPInd) + SPt * 3.*SPC3(SPInd))) * SPiDist
endif

>>> spduparam
SPdm1 = (0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0)))
SPd0 = (0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0))
SPdp1 = (0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0)))
SPdp2 = (SPt * SPt * SPt * 0.166666666666667d0)
if (SPind == 0) then
    SPd = (/SPdm1, SPd0, SPdp1, SPdp2/)
    SPdm1 = 0.
    SPd0 = sum(BCCoefs0(0,:) * SPd)
    SPdp1 = sum(BCCoefs0(1,:) * SPd)
    SPdp2 = sum(BcCoefs0(2,:) * SPd)       
elseif (SPind == NKnot - 1) then
    SPd = (/SPdm1, SPd0, SPdp1, SPdp2/)
    SPdm1 = sum(BCCoefsNm1(0,:) * SPd)
    SPd0 = sum(BCCoefsNm1(1,:) * SPd)
    SPdp1 = 0.
    SPdp2 = 0.  
elseif (SPind == NKnot - 2) then
    SPd = (/SPdm1, SPd0, SPdp1, SPdp2/)
    SPdm1 = sum(BCCoefsNm2(0,:) * SPd)
    SPd0 = sum(BCCoefsNm2(1,:) * SPd)
    SPdp1 = sum(BCCoefsNm2(2,:) * SPd)
    SPdp2 = 0.    
endif
if (CalcDUParam) then
    if (SPdm1 .ne. 0.) DU_Knots(SPInd - 1) = SPdm1 * SPCoef
    if (SPd0  .ne. 0.) DU_Knots(SPInd    ) = SPd0 * SPCoef
    if (SPdp1 .ne. 0.) DU_Knots(SPInd + 1) = SPdp1 * SPCoef
    if (SPdp2 .ne. 0.) DU_Knots(SPInd + 2) = SPdp2 * SPCoef  
endif 

>>> sppairdwparam
if (Thisx >= SPMin .and. Thisx <= SPMax) then
    SPWCoef = DIJ * SPiDist
    SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * SPWCoef
    SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * SPWCoef
    SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * SPWCoef
    SPdp2 = (SPt * SPt * 0.5d0) * SPWCoef
    if (SPind == 0) then
        SPd = (/SPdm1, SPd0, SPdp1, SPdp2/)
        SPdm1 = 0.
        SPd0 = sum(BCCoefs0(0,:) * SPd)
        SPdp1 = sum(BCCoefs0(1,:) * SPd)
        SPdp2 = sum(BcCoefs0(2,:) * SPd)      
    elseif (SPind == NKnot - 1) then
        SPd = (/SPdm1, SPd0, SPdp1, SPdp2/)
        SPdm1 = sum(BCCoefsNm1(0,:) * SPd)
        SPd0 = sum(BCCoefsNm1(1,:) * SPd)
        SPdp1 = 0.
        SPdp2 = 0.
    elseif (SPind == NKnot - 2) then
        SPd = (/SPdm1, SPd0, SPdp1, SPdp2/)
        SPdm1 = sum(BCCoefsNm2(0,:) * SPd)
        SPd0 = sum(BCCoefsNm2(1,:) * SPd)
        SPdp1 = sum(BCCoefsNm2(2,:) * SPd)
        SPdp2 = 0.     
    endif
    if (SPdm1 .ne. 0.) DW_Knots(SPInd - 1) = SPdm1 * SPCoef
    if (SPd0  .ne. 0.) DW_Knots(SPInd    ) = SPd0 * SPCoef
    if (SPdp1 .ne. 0.) DW_Knots(SPInd + 1) = SPdp1 * SPCoef 
    if (SPdp2 .ne. 0.) DW_Knots(SPInd + 2) = SPdp2 * SPCoef
endif    
"""

    SourceCyclic = """
>>> defs
float SPt
int SPInd
int ind1
int ind2
int ind3
int ind4
float SPx
float SPdm1
float SPd0
float SPdp1
float SPdp2
float SPCoef
float Thisx
float SPWCoef

>>> spindex
SPx = (Thisx - SPMin) * SPiDist
SPx = modulo(SPx, real(NInterval))
SPind = int(SPx)
SPt = SPx - real(SPind)

>>> spenergy
THISU = SPCoef * (SPC0(SPInd) + SPt * (SPC1(SPInd) + SPt * (SPC2(SPInd) + SPt * SPC3(SPInd))))

>>> spvirial
THISW = SPCoef * (SPC1(SPInd) + SPt * (2.*SPC2(SPInd) + SPt * 3.*SPC3(SPInd))) * SPx

>>> spdenergy
ThisdU = SPCoef * (SPC1(SPInd) + SPt * (2.*SPC2(SPInd) + SPt * 3.*SPC3(SPInd))) * SPiDist

>>> spduparam
SPdm1 = (0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0)))
SPd0 = (0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0))
SPdp1 = (0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))) 
SPdp2 = (SPt * SPt * SPt * 0.166666666666667d0) 
ind1 = modulo(SPInd - 1, NKnot)
ind2 = modulo(SPInd    , NKnot)
ind3 = modulo(SPInd + 1, NKnot)
ind4 = modulo(SPInd + 2, NKnot)
if (CalcDUParam) then
    DU_Knots(ind1) = SPdm1 * SPCoef
    DU_Knots(ind2) = SPd0 * SPCoef
    DU_Knots(ind3) = SPdp1 * SPCoef
    DU_Knots(ind4) = SPdp2 * SPCoef
endif

>>> sppairdwparam
SPWCoef = DIJ * SPiDist 
SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * SPWCoef
SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * SPWCoef
SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * SPWCoef
SPdp2 = (SPt * SPt * 0.5d0) * SPWCoef
ind1 = modulo(SPInd - 1, NKnot)
ind2 = modulo(SPInd    , NKnot)
ind3 = modulo(SPInd + 1, NKnot)
ind4 = modulo(SPInd + 2, NKnot)
DW_Knots(ind1) = SPdm1 * SPCoef
DW_Knots(ind2) = SPd0 * SPCoef
DW_Knots(ind3) = SPdp1 * SPCoef
DW_Knots(ind4) = SPdp2 * SPCoef
"""    

    
    Names = ["splinebase"]
    #this is a spline
    IsSpline = True
    #can we shift this spline globally in energy without effect?
    KnotsShiftable = False
    #are the knots cyclic?    
    KnotsCyclic = False
    #whether or not to use energy fluctuations (Experimental)
    UsesFluct = False


    def __init__(self, Sys, Label = "", Cut = None, Min = None, Max = None, 
                 Filter = None, Fixed = False, Source = None,
                 InnerBC = None, OuterBC = None, Knots = None, NKnot = 20, 
                 NType = 1, TypeLabels = None):
        #initialize the base class
        PotentialClass.__init__(self, Sys, Label, Cut = Cut, Filter = Filter, 
                                NType = NType, TypeLabels = TypeLabels)
        self.ModuleVars += ["SPDist", "SPiDist", "SPMin", "SPMax", "SPRange", 
                            "SPC0", "SPC1", "SPC2", "SPC3", "SPCoef", "NKnot",
                            "NInterval", "BCCoefs0", "BCCoefsNm1", "BCCoefsNm2"]
        
        #add the knots to the parameters
        if Knots is None:
            self.NKnot = NKnot
            self.Param.Add("Knots", self.NKnot, Value = 0., Fixed = Fixed,
                           Scale = self.Sys.Units.EScale)
        else:
            self.NKnot = len(Knots)
            self.Param.Add("Knots", self.NKnot, Value = Knots, Fixed = Fixed,
                           Scale = self.Sys.Units.EScale)       
        #check that we have enough knots to do interpolation
        if self.NKnot < 4:
            raise ValueError("Must have four or more knot values in potential %s." % self.Label)
            
        #check the boundary conditions for the splines
        #0 = constant value of zero
        #1 = slope of zero
        #2 = zero curvature (default)
        #3 = constant value of zero and slope zero (only for OuterBC)
        if self.KnotsCyclic:
            if not (InnerBC == None and OuterBC == None):
                raise ValueError("Cannot specify boundary conditions for cyclic splines.")
            self.NInterval = self.NKnot
        else:
            if InnerBC not in [0, 1, 2]:
                raise ValueError("Inner boundary condition must be 0, 1, or 2")
            if OuterBC not in [0, 1, 2, 3]:
                raise ValueError("Outer boundary condition must be 0, 1, 2, or 3")                
            if InnerBC is None: InnerBC = 2
            if OuterBC is None: OuterBC = 2
            self.InnerBC = InnerBC
            self.OuterBC = OuterBC
            if OuterBC == 3:
                self.NInterval = self.NKnot
            else:
                self.NInterval = self.NKnot - 1
        
        #make the parameter arrays for the spline
        self.SPCoef = 1.
        if Cut and Min is None and Max is None:
            self.SPMin = 0.
            self.SPMax = Cut
        else:
            self.SPMin = Min
            self.SPMax = Max
        self.SPRange = self.SPMax - self.SPMin
        self.SPDist = 0.
        self.SPiDist = 0.
        self.SPC0 = np.zeros(self.NInterval, dtype=float)
        self.SPC1 = np.zeros(self.NInterval, dtype=float)
        self.SPC2 = np.zeros(self.NInterval, dtype=float)
        self.SPC3 = np.zeros(self.NInterval, dtype=float)
        self.BCCoefs0 = np.zeros((3,4), dtype=float)
        self.BCCoefsNm1 = np.zeros((3,4), dtype=float)
        self.BCCoefsNm2 = np.zeros((3,4), dtype=float)

        #default options for treating spline knots in unsampled regions
        #These are only used by srel optimization
        self.EneInner = None
        self.EneOuter = None
        self.EneInnerInit = None
        self.EneOuterInit = None
        self.EneSlopeInner = None
        self.EneSlopeOuter = None
        self.EneSlopeInnerInit = None
        self.EneSlopeOuterInit = None        
        self.ForceInner = False
        self.ForceOuter = False
        self.KnotMinHistFrac = 0.01
        self.KnotMinHistFracInner = 0
        self.KnotMinHistFracOuter = 0
        self.EneSlopeInnerMin = None
        self.EneSlopeOuterMin = None
        self.EneSlopeBuffer = "0.1kT"
        self.ConstrainSpline = True
        
        #set the source code
        if Source is None: Source = ""
        if self.KnotsCyclic:
            Source = SplinePotentialClass.SourceCyclic + "\n\n" + Source
        else:
            Source = SplinePotentialClass.SourceAcyclic + "\n\n" + Source
        self.SetSource(Source)      


    def Update(self):
        """Updates the potential."""
        n = len(self.Knots)
        #compute the distances
        self.SPDist = self.SPRange / float(self.NInterval)
        self.SPiDist = 1. / self.SPDist     
        #make a shortcut array
        Y = np.zeros(self.NInterval + 3, float)
        #fill in the y-values
        Y[1:1+n] = self.Knots
        if self.KnotsCyclic:
            #make Y vals on boundaries cyclic
            Y[0] = self.Knots[-1]
            Y[-2] = self.Knots[0]
            Y[-1] = self.Knots[1]
        else:
            self.BCCoefs0.fill(0.)
            self.BCCoefsNm1.fill(0.)
            self.BCCoefsNm2.fill(0.)
            #treat the inner boundary condition
            if self.InnerBC == 0:
                #zero value boundary condition: k(-1) = -4 k(0) - k(1)
                self.BCCoefs0[0] = [-4., 1., 0., 0.]
                self.BCCoefs0[1] = [-1., 0., 1., 0.]
                self.BCCoefs0[2] = [0., 0., 0., 1.]
                Y[0] = -4 * Y[1] - Y[2]
            elif self.InnerBC == 1:
                #zero slope boundary condition: k(-1) = k(1)
                self.BCCoefs0[0] = [0., 1., 0., 0.]
                self.BCCoefs0[1] = [1., 0., 1., 0.]
                self.BCCoefs0[2] = [0., 0., 0., 1.]
                Y[0] = Y[2]
            else:
                #zero curvature boundary condition: k(-1) = 2 k(0) - k(1)
                self.BCCoefs0[0] = [2., 1., 0., 0.]
                self.BCCoefs0[1] = [-1., 0., 1., 0.]
                self.BCCoefs0[2] = [0., 0., 0., 1.]
                Y[0] = 2 * Y[1] - Y[2]
            if self.OuterBC == 0:
                #zero value boundary condition: k(N) = -k(N-2) - 4 k(N-1)
                self.BCCoefsNm2[0] = [1., 0., 0., 0.]
                self.BCCoefsNm2[1] = [0., 1., 0., -1.]
                self.BCCoefsNm2[2] = [0., 0., 1., -4.]
                self.BCCoefsNm1[0] = [1., 0., -1., 0.]
                self.BCCoefsNm1[1] = [0., 1., -4., 0.]
                Y[n+1] = -Y[n-1] - 4 * Y[n]
            elif self.OuterBC == 1:
                #zero slope boundary condition: k(N) = k(N-2)
                self.BCCoefsNm2[0] = [1., 0., 0., 0.]
                self.BCCoefsNm2[1] = [0., 1., 0., 1.]
                self.BCCoefsNm2[2] = [0., 0., 1., 0.]
                self.BCCoefsNm1[0] = [1., 0., 1., 0.]
                self.BCCoefsNm1[1] = [0., 1., 0., 0.]
                Y[n+1] = Y[n-1]
            elif self.OuterBC == 3:
                #zero value and slope boundary condition: k(N) = -0.5 k(N-1) and k(N+1) = k(N-1)
                self.BCCoefsNm2[0] = [1., 0., 0., 0.]
                self.BCCoefsNm2[1] = [0., 1., 0., 0.]
                self.BCCoefsNm2[2] = [0., 0., 1., -0.5]
                self.BCCoefsNm1[0] = [1., 0., 0., 0.]
                self.BCCoefsNm1[1] = [0., 1., -0.5, 1.]
                Y[n+1] = -0.5 * Y[n]
                Y[n+2] = Y[n]
            else:
                #zero curvature boundary condition: k(N) = 2 k(N-1) - k(N-2)
                self.BCCoefsNm2[0] = [1., 0., 0., 0.]
                self.BCCoefsNm2[1] = [0., 1., 0., -1.]
                self.BCCoefsNm2[2] = [0., 0., 1., 2.]
                self.BCCoefsNm1[0] = [1., 0., -1., 0.]
                self.BCCoefsNm1[1] = [0., 1., 2., 0.]   
                Y[n+1] = 2 * Y[n] - Y[n-1]
        #make the coefficients
        C = np.zeros((4, self.NInterval), float)
        for i in range(self.NInterval):
            C[:,i] = np.dot(self.BFunc, Y[i:i+4])
        self.SPC0 = C[0]
        self.SPC1 = C[1]
        self.SPC2 = C[2]
        self.SPC3 = C[3] 
        
        
    def RelaxKnotConstraints(self):
        """Sets EneInner, EneSlopeInner, EneOuter, and EneSlopeOuter to None."""
        #These are only used by srel optimization
        self.EneInner = None
        self.EneOuter = None
        self.EneInnerInit = None
        self.EneOuterInit = None
        self.EneSlopeInner = None
        self.EneSlopeOuter = None
        self.ForceInner = False
        self.ForceOuter = False

        
    def SetParam(self, Knots):
        """Sets parameters for this potential."""
        try:
            self.Knots = Knots
        except ValueError:
            raise ValueError("Improper value for Knots. Must be single float or length-%d" % self.NKnot)
        self.Update()        


    def FitSpline(self, ArgVals, EneVals, Weights = None, MaxKnot = None, MinKnot = None):
        """Fits the knot points using least squares regression.
ArgVals: values of arguments to the energy function (eg, pair distances0
EneVals: target energies at these argument values  
Weights: weights for different parts of the fit (optional)"""
        from sim.mathutil import LeastSquaresFit
        if not len(ArgVals) == len(EneVals):
            raise ValueError("ArgVals and EneVals are different lengths.")
        if not Weights is None:
            if not len(ArgVals) == len(Weights):
                raise ValueError("ArgVals and Weights are different lengths.")
        #make inputs to least squares regression
        Xs = []
        for arg in ArgVals:
            x = self.DUParam(arg)
            if not np.all(x==0.):
                Xs.append(x)
        Xs = np.array(Xs, dtype=float)
        Ys = np.array(EneVals, dtype=float)   
        #now estimate the spline knots
        Mask = np.any(Xs != 0., axis = 0)
        self.Knots[:] = 0.
        self.Knots[Mask] = LeastSquaresFit(Xs[:,Mask], Ys, Weights)
        #filter for maximum energy
        if not MaxKnot is None:
            for (i,v) in enumerate(self.Knots):
                if v > MaxKnot:
                    self.Knots[i] = MaxKnot 
        #filter for minimum energy
        if not MinKnot is None:
            for (i,v) in enumerate(self.Knots):
                if v < MinKnot:
                    self.Knots[i] = MinKnot                     
        self.Update()
              
      
    def Estimate(self, MinKnot = None, MaxKnot = None):
        """Estimates spline potential parameters based on arg counts and histogram."""           
        kT = self.Sys.TempSet * self.Sys.Units.kB
        ArgVals = np.zeros(self.Arg.HistNBin, float)
        EneVals = np.zeros(self.Arg.HistNBin, float)
        Weights = np.zeros(self.Arg.HistNBin, float)
        HistSum = np.sum(self.Arg.Hist[0])
        if HistSum == 0:
            raise ValueError("No entries in histogram.")
        HistMin = np.min(self.Arg.Hist[0][self.Arg.Hist[0] > 0]) 
        HistErr = self.HistFracErr * HistMin / HistSum
        for (i, h) in enumerate(self.Arg.Hist[0]):
            x = self.Arg.HistMin[0] + self.Arg.HistBinw[0] * (i + 0.5)
            ArgVals[i] = x
            w = HistErr + h / HistSum
            EneVals[i] = -kT * np.log(w)
            Weights[i] = w
        self.FitSpline(ArgVals, EneVals, Weights, MaxKnot, MinKnot)


    def KnotHist(self):
        """Returns a histogram of frequencies in which different knot values are used."""
        KHist = np.zeros(len(self.Knots), float)
        for (i, h) in enumerate(self.Arg.Hist[0]):
            x = self.Arg.HistMin[0] + self.Arg.HistBinw[0] * (i + 0.5)
            KHist += self.DUParam(x) * h
        return KHist
    
    
    def KnotInd(self, x):
        """Returns the knot index for this potential."""
        if self.KnotsCyclic:
            x = np.mod(x - self.SPMin, self.SPRange) * self.SPiDist
        else:
            if x >= self.SPMax or x < self.SPMin:
                return None
            x = (x - self.SPMin) * self.SPiDist
        return min(int(x), self.NInterval - 1)    
        
    
    def Val(self, x):
        """Returns the value of the potential."""
        if self.KnotsCyclic:
            x = np.mod(x - self.SPMin, self.SPRange) * self.SPiDist
        else:
            x = (x - self.SPMin) * self.SPiDist
        i = max(min(int(x), self.NInterval - 1), 0)
        t = max(min(x - float(i), 1.0), 0.0)
        return (self.SPC0[i] + t * (self.SPC1[i] + t * (self.SPC2[i] + t * self.SPC3[i]))) * self.SPCoef
        
    
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if self.KnotsCyclic:
            x = np.mod(x - self.SPMin, self.SPRange) * self.SPiDist
        else:
            x = (x - self.SPMin) * self.SPiDist
        if x < 0. or x >= self.NInterval:
            return 0.
        i = max(min(int(x), self.NInterval - 1), 0)
        t = max(min(x - float(i), 1.0), 0.0)
        return (self.SPC1[i] + t * (2.*self.SPC2[i] + t * 3.*self.SPC3[i])) * self.SPiDist * self.SPCoef


    def DUParam(self, x):
        """Returns the values of the derivative with respect to spline knots."""
        if self.KnotsCyclic:
            x = np.mod(x - self.SPMin, self.SPRange) * self.SPiDist
        else:
            x = (x - self.SPMin) * self.SPiDist
        SPInd = max(min(int(x), self.NInterval - 1), 0)
        t = max(min(x - float(SPInd), 1.0), 0.0)
        d = np.zeros(self.NKnot, float)
        SPdm1 = 0.166666666666667 + t * (-0.5 + t * (0.5 - t * 0.166666666666667))
        SPd0 = 0.666666666666667 + t * t * (-1.0 + t * 0.5)
        SPdp1 = 0.166666666666667 + t * (0.5 + t * (0.5 - t * 0.5))
        SPdp2 = t * t * t * 0.166666666666667
        if self.KnotsCyclic:
            if SPInd == 0:
                d[self.NKnot - 1] = SPdm1
                d[SPInd    ] = SPd0
                d[SPInd + 1] = SPdp1
                d[SPInd + 2] = SPdp2
            elif SPInd == self.NKnot - 1:
                d[SPInd - 1] = SPdm1
                d[SPInd    ] = SPd0
                d[0        ] = SPdp1
                d[1        ] = SPdp2
            elif SPInd == self.NKnot - 2:
                d[SPInd - 1] = SPdm1
                d[SPInd    ] = SPd0
                d[SPInd + 1] = SPdp1
                d[0        ] = SPdp2
            else:
                d[SPInd - 1] = SPdm1
                d[SPInd    ] = SPd0
                d[SPInd + 1] = SPdp1
                d[SPInd + 2] = SPdp2
        else:
            SPd = np.array([SPdm1, SPd0, SPdp1, SPdp2])
            if SPInd == 0:
                d[SPInd    ] = np.sum(self.BCCoefs0[0,:] * SPd)
                d[SPInd + 1] = np.sum(self.BCCoefs0[1,:] * SPd)
                d[SPInd + 2] = np.sum(self.BCCoefs0[2,:] * SPd)
            elif SPInd == self.NKnot - 1:
                d[SPInd - 1] = np.sum(self.BCCoefsNm1[0,:] * SPd)
                d[SPInd    ] = np.sum(self.BCCoefsNm1[1,:] * SPd)
            elif SPInd == self.NKnot - 2:
                d[SPInd - 1] = np.sum(self.BCCoefsNm2[0,:] * SPd)
                d[SPInd    ] = np.sum(self.BCCoefsNm2[1,:] * SPd)
                d[SPInd + 1] = np.sum(self.BCCoefsNm2[2,:] * SPd)
            else:
                d[SPInd - 1] = SPdm1
                d[SPInd    ] = SPd0
                d[SPInd + 1] = SPdp1
                d[SPInd + 2] = SPdp2
        d = d * self.SPCoef
        return d

   