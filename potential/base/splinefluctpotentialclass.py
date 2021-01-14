### Class definitions for spline potential energy functions in SIM suite.
### coded by MSS


import numpy as np
from splinepotentialclass import SplinePotentialClass


class SplineFluctPotentialClass(SplinePotentialClass):
    """Generic class for developing spline-based potential energy functions with energy fluctuations."""
    
    #bspline matrix
    BFunc = np.array([[ 1.,  4.,  1.,  0.],
                      [-3.,  0.,  3.,  0.],
                      [ 3., -6.,  3.,  0.],
                      [-1.,  3., -3.,  1.]], float) / 6.
    
    
    SourceCommon = """
>>> defs
float FLt
int FLInd
float FLx
float FLdm1
float FLd0
float FLdp1
float FLdp2
float Thisx
float FLd(4)
float iFluctA0
float DeltaE
float val1
float val2
float val3
float val4
float val5
float FluctkT

>>> fluct
ThisA = FLC0(FLInd) + FLt * (FLC1(FLInd) + FLt * (FLC2(FLInd) + FLt * FLC3(FLInd)))
FluctA0 = FluctA0 + ThisA

>>> init
if (CalcFluct) then
    FLg = 0.
    FLf = 0.
endif

>>> final
if (CalcFluct) then
    !calculate contributions to energy fluctuations
    FluctkT = 1. / FluctBeta
    FluctE0 = PEnergy - 0.5 * FluctBeta * FluctA0
    !calculate shortcut variables
    iFluctA0 = 1. / FluctA0
    DeltaE = FluctE - FluctE0
    !calculate -kB T ln P(E)
    FluctTerm = FluctkT * (0.5 * DeltaE**2 * iFluctA0 + 0.5 * log(2.*pi*FluctA0))
    val1 = -FluctkT * DeltaE * iFluctA0 
    val2 = FluctkT * (FluctA0 - DeltaE**2) * 0.5 * iFluctA0**2 + DeltaE * 0.5 * iFluctA0
    val3 = FluctkT * iFluctA0
    val4 = -0.5 * iFluctA0 + FluctkT * DeltaE * iFluctA0**2
    val5 = -FluctkT * (FluctA0 - 2. * DeltaE**2) * 0.5 * iFluctA0**3 - DeltaE * iFluctA0**2 + 0.25 * FluctBeta * iFluctA0
    !first derivative of LnP with respect to interaction knots
    do i = 0, NKnot - 1
        DU_Knots(i) = val1 * FLf(i)
    enddo
    !first derivative of LnP with respect to fluct knots
    do i = 0, FLNKnot - 1
        DU_FLKnots(i) = val2 * FLg(i)
    enddo
    !second derivative of LnP with respect to interactions
    do i = 0, NKnot - 1
        do j = i, NKnot - 1
            DDU_Knots_Knots(i,j) = val3 * FLf(i) * FLf(j)
        enddo
    enddo
    !second cross derivative with respect to interactions and fluct
    do i = 0, NKnot - 1
        do j = 0, FLNKnot - 1
            DDU_Knots_FLKnots(i,j) = val4 * FLf(i) * FLg(j)
        enddo
    enddo
    !second derivative term with respect to fluct
    do i = 0, FLNKnot - 1
        do j = i, FLNKnot - 1
            DDU_FLKnots_FLKnots(i,j) = val5 * FLg(i) * FLf(j)
        enddo
    enddo
endif
"""

    SourceAcyclic = SourceCommon + """
>>> fluctindex
FLx = (Thisx - FLMin) * FLiDist
FLind = max(min(int(FLx), FLNInterval - 1), 0) 
FLt = max(min(FLx - real(FLind), 1.d0), 0.d0)

>>> spduparam
if (CalcFluct) then 
    !store these values to calculate fluctuation term to derivatives (at end)
    if (SPdm1 .ne. 0) FLf(SPInd - 1) = FLf(SPInd - 1) + SPdm1 
    if (SPd0  .ne. 0) Flf(SPInd    ) = FLf(SPInd) + SPd0 
    if (SPdp1 .ne. 0) Flf(SPInd + 1) = FLf(SPInd + 1) + SPdp1 
    if (SPdp2 .ne. 0) Flf(SPInd + 2) = FLf(SPInd + 2) + SPdp2 
endif

>>> fluctdparam
FLdm1 = (0.166666666666667d0 + FLt * (-0.5d0 + FLt * (0.5d0 - FLt * 0.166666666666667d0)))
FLd0 = (0.666666666666667d0 + FLt * FLt * (-1.0d0 + FLt * 0.5d0))
FLdp1 = (0.166666666666667d0 + FLt * (0.5d0 + FLt * (0.5d0 - FLt * 0.5d0)))
FLdp2 = (FLt * FLt * FLt * 0.166666666666667d0)
if (FLind == 0) then
    FLd = (/FLdm1, FLd0, FLdp1, FLdp2/)
    FLdm1 = 0.
    FLd0 = sum(FLBCCoefs0(0,:) * FLd)
    FLdp1 = sum(FLBCCoefs0(1,:) * FLd)
    FLdp2 = sum(FLBCCoefs0(2,:) * FLd)    
elseif (FLind == FLNKnot - 1) then
    FLd = (/FLdm1, FLd0, FLdp1, FLdp2/)
    FLdm1 = sum(FLBCCoefsNm1(0,:) * FLd)
    FLd0 = sum(FLBCCoefsNm1(1,:) * FLd)
    FLdp1 = 0.
    FLdp2 = 0.   
elseif (FLind == FLNKnot - 2) then
    FLd = (/FLdm1, FLd0, FLdp1, FLdp2/)
    FLdm1 = sum(FLBCCoefsNm2(0,:) * FLd)
    FLd0 = sum(FLBCCoefsNm2(1,:) * FLd)
    FLdp1 = sum(FLBCCoefsNm2(2,:) * FLd)
    FLdp2 = 0.     
endif
if (FLdm1 .ne. 0.) FLg(FLInd - 1) = FLg(FLInd - 1) + FLdm1 
if (FLd0  .ne. 0.) FLg(FLInd    ) = FLg(FLInd    ) + FLd0 
if (FLdp1 .ne. 0.) FLg(FLInd + 1) = FLg(FLInd + 1) + FLdp1   
if (FLdp2 .ne. 0.) FLg(FLInd + 2) = FLg(FLInd + 2) + FLdp2 
"""

    SourceCyclic = SourceCommon + """
>>> defs
int ind1
int ind2
int ind3
int ind4

>>> fluctindex
FLx = (Thisx - FLMin) * FLiDist
FLx = modulo(FLx, real(FLNInterval))
FLind = int(FLx)
FLt = FLx - real(FLind)

>>> spduparam
if (CalcFluct) then
    !store these values to calculate fluctuation term to derivatives (at end)
    FLf(ind1) = SPdm1 * SPCoef
    Flf(ind2) = SPd0 * SPCoef
    Flf(ind3) = SPdp1 * SPCoef
    Flf(ind4) = SPdp2 * SPCoef
endif

>>> fluctdparam
FLdm1 = (0.166666666666667d0 + FLt * (-0.5d0 + FLt * (0.5d0 - FLt * 0.166666666666667d0)))
FLd0 = (0.666666666666667d0 + FLt * FLt * (-1.0d0 + FLt * 0.5d0))
FLdp1 = (0.166666666666667d0 + FLt * (0.5d0 + FLt * (0.5d0 - FLt * 0.5d0))) 
FLdp2 = (FLt * FLt * FLt * 0.166666666666667d0) 
ind1 = modulo(FLInd - 1, NKnot)
ind2 = modulo(FLInd    , NKnot)
ind3 = modulo(FLInd + 1, NKnot)
ind4 = modulo(FLInd + 2, NKnot)
FLg(ind1) = FLg(ind1) + FLdm1
FLg(ind2) = FLg(ind2) + FLd0
FLg(ind3) = FLg(ind3) + FLdp1
FLg(ind4) = FLg(ind4) + FLdp2
"""    

    
    Names = ["splinefluctbase"]
    #whether or not to use energy fluctuations (Experimental)
    UsesFluct = True


    def __init__(self, Sys, Cut = None, FLMin = None, FLMax = None, 
                 FluctMin = None, FluctSoftMin = None,
                 FLInnerBC = 0, FLOuterBC = 0, FLNKnot = 5, FLKnots = None, 
                 Source = None, **kwargs):
        #initialize the base class
        SplinePotentialClass.__init__(self, Sys, Cut = Cut, **kwargs)
        
        self.ModuleVars += ["FLDist", "FLiDist", "FLMin", "FLMax", 
                            "FLC0", "FLC1", "FLC2", "FLC3", "FLf", "FLg",
                            "FLNKnot", "FLNInterval",
                            "FLBCCoefs0", "FLBCCoefsNm1", "FLBCCoefsNm2"]
        
        #check scales and minimum values
        FLScale = self.Sys.Units.EScale**2
        #minimum value of the fluctuations at any distance
        if FluctMin is None:
            FluctMin = 0.01 * FLScale
        if FluctSoftMin is None:
            FluctSoftMin = 2 * FluctMin
        self.FluctMin = FluctMin
        self.FluctSoftMin = FluctSoftMin
        #add the knot parameters
        if FLKnots is None:
            self.FLNKnot = FLNKnot
            self.Param.Add("FLKnots", self.FLNKnot, Value = FLScale, Scale = FLScale,
                           SoftMin = FluctSoftMin, Min = FluctMin)
        else:
            self.FLNKnot = len(FLKnots)
            self.Param.Add("FLKnots", self.FLNKnot, Value = FLKnots, Scale = FLScale,
                           SoftMin = FluctSoftMin, Min = FluctMin)
        #check that we have enough knots to do interpolation
        if self.FLNKnot < 4:
            raise ValueError("Must have four or more knot values in fluct potential %s." % self.Label)
            
        #check the boundary conditions for the splines
        #0 = constant value of zero
        #1 = slope of zero
        #2 = zero curvature (default)
        #3 = constant value of zero and slope zero (only for FLOuterBC)
        if self.KnotsCyclic:
            if not (FLInnerBC == None and FLOuterBC == None):
                raise ValueError("Cannot specify boundary conditions for cyclic splines.")
            self.FLNInterval = self.FLNKnot
        else:
            if FLInnerBC not in [0, 1, 2]:
                raise ValueError("Inner boundary condition must be 0, 1, or 2")
            if FLOuterBC not in [0, 1, 2, 3]:
                raise ValueError("Outer boundary condition must be 0, 1, 2, or 3")                
            if FLInnerBC is None: FLInnerBC = 2
            if FLOuterBC is None: FLOuterBC = 2
            self.FLInnerBC = FLInnerBC
            self.FLOuterBC = FLOuterBC
            if FLOuterBC == 3:
                self.FLNInterval = self.FLNKnot
            else:
                self.FLNInterval = self.FLNKnot - 1
        
        #make the parameter arrays for the spline
        if Cut and FLMin is None and FLMax is None:
            self.FLMin = 0.
            self.FLMax = Cut
        else:
            self.FLMin = FLMin
            self.FLMax = FLMax
        self.FLRange = self.FLMax - self.FLMin
        self.FLDist = 0.
        self.FLiDist = 0.
        self.FLC0 = np.zeros(self.FLNInterval, dtype=float)
        self.FLC1 = np.zeros(self.FLNInterval, dtype=float)
        self.FLC2 = np.zeros(self.FLNInterval, dtype=float)
        self.FLC3 = np.zeros(self.FLNInterval, dtype=float)
        self.FLBCCoefs0 = np.zeros((3,4), dtype=float)
        self.FLBCCoefsNm1 = np.zeros((3,4), dtype=float)
        self.FLBCCoefsNm2 = np.zeros((3,4), dtype=float)
        #arrays for counting uses of knot functions
        self.FLf = np.zeros(self.NKnot, dtype=float)
        self.FLg = np.zeros(self.FLNKnot, dtype=float)

        #default options for treating spline knots in unsampled regions
        #These are only used by srel optimization
        #0 = make zero in all unsampled regions
        #1 = make constant slope in unsampled regions
        self.FluctTreatment = 0  
        
        #get the source code
        if Source is None: Source = ""
        if self.KnotsCyclic:
            Source = SplineFluctPotentialClass.SourceCyclic + "\n\n" + Source
        else:
            Source = SplineFluctPotentialClass.SourceAcyclic + "\n\n" + Source
        self.SetSource(Source)     


    def Update(self):
        """Updates the potential."""
        SplinePotentialClass.Update(self)
        n = len(self.FLKnots)
        #compute the distances
        self.FLDist = self.FLRange / float(self.FLNInterval)
        self.FLiDist = 1. / self.FLDist     
        #make a shortcut array
        Y = np.zeros(self.FLNInterval + 3, float)
        #fill in the y-values
        Y[1:1+n] = self.FLKnots
        if self.KnotsCyclic:
            #make Y vals on boundaries cyclic
            Y[0] = self.FLKnots[-1]
            Y[-2] = self.FLKnots[0]
            Y[-1] = self.FLKnots[1]
        else:
            self.FLBCCoefs0.fill(0.)
            self.FLBCCoefsNm1.fill(0.)
            self.FLBCCoefsNm2.fill(0.)
            #treat the inner boundary condition
            if self.FLInnerBC == 0:
                #zero value boundary condition: k(-1) = -4 k(0) - k(1)
                self.FLBCCoefs0[0] = [-4., 1., 0., 0.]
                self.FLBCCoefs0[1] = [-1., 0., 1., 0.]
                self.FLBCCoefs0[2] = [0., 0., 0., 1.]
                Y[0] = -4 * Y[1] - Y[2]
            elif self.FLInnerBC == 1:
                #zero slope boundary condition: k(-1) = k(1)
                self.FLBCCoefs0[0] = [0., 1., 0., 0.]
                self.FLBCCoefs0[1] = [1., 0., 1., 0.]
                self.FLBCCoefs0[2] = [0., 0., 0., 1.]
                Y[0] = Y[2]
            else:
                #zero curvature boundary condition: k(-1) = 2 k(0) - k(1)
                self.FLBCCoefs0[0] = [2., 1., 0., 0.]
                self.FLBCCoefs0[1] = [-1., 0., 1., 0.]
                self.FLBCCoefs0[2] = [0., 0., 0., 1.]
                Y[0] = 2 * Y[1] - Y[2]
            if self.FLOuterBC == 0:
                #zero value boundary condition: k(N) = -k(N-2) - 4 k(N-1)
                self.FLBCCoefsNm2[0] = [1., 0., 0., 0.]
                self.FLBCCoefsNm2[1] = [0., 1., 0., -1.]
                self.FLBCCoefsNm2[2] = [0., 0., 1., -4.]
                self.FLBCCoefsNm1[0] = [1., 0., -1., 0.]
                self.FLBCCoefsNm1[1] = [0., 1., -4., 0.]
                Y[n+1] = -Y[n-1] - 4 * Y[n]
            elif self.FLOuterBC == 1:
                #zero slope boundary condition: k(N) = k(N-2)
                self.FLBCCoefsNm2[0] = [1., 0., 0., 0.]
                self.FLBCCoefsNm2[1] = [0., 1., 0., 1.]
                self.FLBCCoefsNm2[2] = [0., 0., 1., 0.]
                self.FLBCCoefsNm1[0] = [1., 0., 1., 0.]
                self.FLBCCoefsNm1[1] = [0., 1., 0., 0.]
                Y[n+1] = Y[n-1]
            elif self.FLOuterBC == 3:
                #zero value and slope boundary condition: k(N) = -0.5 k(N-1) and k(N+1) = k(N-1)
                self.FLBCCoefsNm2[0] = [1., 0., 0., 0.]
                self.FLBCCoefsNm2[1] = [0., 1., 0., 0.]
                self.FLBCCoefsNm2[2] = [0., 0., 1., -0.5]
                self.FLBCCoefsNm1[0] = [1., 0., 0., 0.]
                self.FLBCCoefsNm1[1] = [0., 1., -0.5, 1.]
                Y[n+1] = -0.5 * Y[n]
                Y[n+2] = Y[n]
            else:
                #zero curvature boundary condition: k(N) = 2 k(N-1) - k(N-2)
                self.FLBCCoefsNm2[0] = [1., 0., 0., 0.]
                self.FLBCCoefsNm2[1] = [0., 1., 0., -1.]
                self.FLBCCoefsNm2[2] = [0., 0., 1., 2.]
                self.FLBCCoefsNm1[0] = [1., 0., -1., 0.]
                self.FLBCCoefsNm1[1] = [0., 1., 2., 0.]   
                Y[n+1] = 2 * Y[n] - Y[n-1]
        #make the coefficients
        C = np.zeros((4, self.FLNInterval), float)
        for i in range(self.FLNInterval):
            C[:,i] = np.dot(self.BFunc, Y[i:i+4])
        self.FLC0 = C[0]
        self.FLC1 = C[1]
        self.FLC2 = C[2]
        self.FLC3 = C[3] 
        
                
    def SetParam(self, Knots, FLKnots):
        """Sets parameters for this potential."""
        SplinePotentialClass.SetParam(Knots)
        try:
            self.FLKnots = FLKnots
        except ValueError:
            raise ValueError("Improper value for FLKnots. Must be single float or length-%d" % self.FLNKnot)
        self.Update()        


    def FLKnotHist(self):
        """Returns a histogram of frequencies in which different fluctuation knot values are used."""
        KHist = np.zeros(len(self.FLKnots), float)
        for (i, h) in enumerate(self.Arg.Hist[0]):
            x = self.Arg.HistMin[0] + self.Arg.HistBinw[0] * (i + 0.5)
            KnotInd = self.FLKnotInd(x)
            if KnotInd is None:
                continue
            KHist[KnotInd] += h
            if KnotInd + 1 < self.FLNKnot:
                KHist[KnotInd+1] += h
        return KHist
    
    
    def FLKnotInd(self, x):
        """Returns the fluctuation knot index for this potential."""
        if self.KnotsCyclic:
            x = np.mod(x - self.FLMin, self.FLRange) * self.FLiDist
        else:
            if x >= self.FLMax or x < self.FLMin:
                return None
            x = (x - self.FLMin) * self.FLiDist
        return min(int(x), self.FLNInterval - 1)    
        
    
    def FLVal(self, x):
        """Returns the value of the fluctuations."""
        if self.KnotsCyclic:
            x = np.mod(x - self.FLMin, self.FLRange) * self.FLiDist
        else:
            x = (x - self.FLMin) * self.FLiDist
        i = max(min(int(x), self.FLNInterval - 1), 0)
        t = max(min(x - float(i), 1.0), 0.0)
        return self.FLC0[i] + t * (self.FLC1[i] + t * (self.FLC2[i] + t * self.FLC3[i]))
        

    def FLDParam(self, x):
        """Returns the values of the derivative with respect to fluctuation spline knots."""
        if self.KnotsCyclic:
            x = np.mod(x - self.FLMin, self.FLRange) * self.FLiDist
        else:
            x = (x - self.FLMin) * self.FLiDist
        SPInd = max(min(int(x), self.FLNInterval - 1), 0)
        t = max(min(x - float(SPInd), 1.0), 0.0)
        d = np.zeros(self.FLNKnot, float)
        FLdm1 = 0.166666666666667 + t * (-0.5 + t * (0.5 - t * 0.166666666666667))
        FLd0 = 0.666666666666667 + t * t * (-1.0 + t * 0.5)
        FLdp1 = 0.166666666666667 + t * (0.5 + t * (0.5 - t * 0.5))
        FLdp2 = t * t * t * 0.166666666666667
        if self.KnotsCyclic:
            if SPInd == 0:
                d[self.FLNKnot - 1] = FLdm1
                d[SPInd    ] = FLd0
                d[SPInd + 1] = FLdp1
                d[SPInd + 2] = FLdp2
            elif SPInd == self.FLNKnot - 1:
                d[SPInd - 1] = FLdm1
                d[SPInd    ] = FLd0
                d[0        ] = FLdp1
                d[1        ] = FLdp2
            elif SPInd == self.FLNKnot - 2:
                d[SPInd - 1] = FLdm1
                d[SPInd    ] = FLd0
                d[SPInd + 1] = FLdp1
                d[0        ] = FLdp2
            else:
                d[SPInd - 1] = FLdm1
                d[SPInd    ] = FLd0
                d[SPInd + 1] = FLdp1
                d[SPInd + 2] = FLdp2
        else:
            FLd = np.array([FLdm1, FLd0, FLdp1, FLdp2])
            if SPInd == 0:
                d[SPInd    ] = np.sum(self.FLBCCoefs0[0,:] * FLd)
                d[SPInd + 1] = np.sum(self.FLBCCoefs0[1,:] * FLd)
                d[SPInd + 2] = np.sum(self.FLBCCoefs0[2,:] * FLd)
            elif SPInd == self.FLNKnot - 1:
                d[SPInd - 1] = np.sum(self.FLBCCoefsNm1[0,:] * FLd)
                d[SPInd    ] = np.sum(self.FLBCCoefsNm1[1,:] * FLd)
            elif SPInd == self.FLNKnot - 2:
                d[SPInd - 1] = np.sum(self.FLBCCoefsNm2[0,:] * FLd)
                d[SPInd    ] = np.sum(self.FLBCCoefsNm2[1,:] * FLd)
                d[SPInd + 1] = np.sum(self.FLBCCoefsNm2[2,:] * FLd)
            else:
                d[SPInd - 1] = FLdm1
                d[SPInd    ] = FLd0
                d[SPInd + 1] = FLdp1
                d[SPInd + 2] = FLdp2
        return d

   