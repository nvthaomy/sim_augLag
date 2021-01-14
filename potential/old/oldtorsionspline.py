#/usr/bin/env python


import numpy as np
import base



class OldTorsionSpline(base.SplinePotentialClass):
    """Spline torsion angle interactions."""

    Names = ["oldtorsionspline"]    

    Source = """
>>> defs
float SPt
int SPInd
float SPx
float SPdm1
float SPd0
float SPdp1
float SPdp2

>>> torsionenergydparam
SPx = (PHI + PI) * SPiDist
SPind = max(min(int(SPx), [[NKnot]] - 1), 0)
SPt = SPx - real(SPind)
ThisU = SPC0(SPInd) + SPt * (SPC1(SPInd) + SPt * (SPC2(SPInd) + SPt * SPC3(SPInd)))
[UPDATEENERGY]
if (CalcDUParam) then
    SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
    SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
    SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
    SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
    if (SPind == 0) then
        DU_Knots([[NKnot]] - 1) = SPdm1
        DU_Knots(0            ) = SPd0
        DU_Knots(1            ) = SPdp1
        DU_Knots(2            ) = SPdp2
    elseif (SPind == [[NKnot]] - 1) then
        DU_Knots([[NKnot]] - 2) = SPdm1
        DU_Knots([[NKnot]] - 1) = SPd0
        DU_Knots(0            ) = SPdp1
        DU_Knots(1            ) = SPdp2
    elseif (SPind == [[NKnot]] - 2) then
        DU_Knots([[NKnot]] - 3) = SPdm1
        DU_Knots([[NKnot]] - 2) = SPd0
        DU_Knots([[NKnot]] - 1) = SPdp1
        DU_Knots(0            ) = SPdp2
    else
        DU_Knots(SPInd - 1) = SPdm1
        DU_Knots(SPInd    ) = SPd0
        DU_Knots(SPInd + 1) = SPdp1
        DU_Knots(SPInd + 2) = SPdp2
    endif
endif

>>> argtorsionloop
ARGVAL = Phi
ARGTYPE = 0
[ARGGET]

>>> argeval
Phi = ARGVAL
[torsionenergydparam]

>>> torsionloop
if (TorDefined) then
    [torsionenergydparam]
    if (CalcForce) then
        ThisDU =  (SPC1(SPInd) + SPt * (2.*SPC2(SPInd) + SPt * 3.*SPC3(SPInd))) * SPiDist 
        Forcei =  (djk * ThisdU) * crossijk / crosssqijk
        Forcel = -(djk * ThisdU) * crossjkl / crosssqjkl
        Forcej = (dot_product(rkl,rjk) / djk**2) * Forcel &
             & - (dot_product(rij,rjk) / djk**2 + 1.d0) * Forcei
        Forcek = (dot_product(rij,rjk) / djk**2) * Forcei &
             & - (dot_product(rkl,rjk) / djk**2 + 1.d0) * Forcel	
        Force(i,:) = Force(i,:) + Forcei
        Force(j,:) = Force(j,:) + Forcej
        Force(k,:) = Force(k,:) + Forcek
        Force(l,:) = Force(l,:) + Forcel
    endif
endif
"""
    
    TestArgs = {"Knots" : [-5, -3, 0, 1, -2]}
    
    Type = base.ptypes.TorsionPotential
    KnotsShiftable = True
    KnotsCyclic = True
    

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False, Knots = None, NKnot = 20):
        """Initializes a spline torsion angle interaction."""
        base.SplinePotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                           Filter = Filter, Fixed = Fixed,
                                           Knots = Knots, NKnot = NKnot,
                                           Source = OldTorsionSpline.Source,
                                           Min = -np.pi, Max = np.pi)            
        self.Update()    

    def SetParam(self, Knots):
        """Sets parameters for this potential."""
        if not len(Knots) == self.NKnot:
            raise IndexError("Require a knot array with %d values." % self.NKnot)
        self.Knots = Knots
        self.Update()

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.Knots.MaxChange = MaxFracChange * self.Sys.Units.EScale
        
    def KnotInd(self, x):
        """Returns the lower knot index for this potential."""
        x = np.mod(x - np.pi, 2*np.pi)
        return min(int(x * self.SPiDist), self.NKnot - 1)    

    def Val(self, x):
        """Returns the value of the potential."""
        x = np.mod(x - np.pi, 2*np.pi)
        x = x * self.SPiDist
        i = min(int(x), self.NKnot - 1)
        t = x - float(i)
        return self.SPC0[i] + t * (self.SPC1[i] + t * (self.SPC2[i] + t * self.SPC3[i]))
        
    def DVal(self, x):
        """Returns the derivative of the potential."""
        x = np.mod(x - np.pi, 2*np.pi)
        x = x * self.SPiDist
        i = min(int(x), self.NKnot - 1)
        t = x - float(i)
        return (self.SPC1[i] + t * (2.*self.SPC2[i] + t * 3.*self.SPC3[i])) * self.SPiDist  
        
    def DUParam(self, x):
        """Returns the values of the derivative with respect to the spline knots."""
        x = np.mod(x - np.pi, 2*np.pi)
        x = x * self.SPiDist
        i = min(int(x), self.NKnot - 1)
        t = x - float(i)
        d = np.zeros(self.NKnot, float)
        SPdm1 = 0.166666666666667 + t * (-0.5 + t * (0.5 - t * 0.166666666666667))
        SPd0 = 0.666666666666667 + t * t * (-1.0 + t * 0.5)
        SPdp1 = 0.166666666666667 + t * (0.5 + t * (0.5 - t * 0.5))
        SPdp2 = t * t * t * 0.166666666666667
        d[(i - 1) % self.NKnot] = SPdm1
        d[ i      % self.NKnot] = SPd0 
        d[(i + 1) % self.NKnot] = SPdp1
        d[(i + 2) % self.NKnot] = SPdp2
        return d       

    def Update(self):
        """Updates the potential."""
        #compute the distances
        self.SPDist = 2.*np.pi / float(self.NKnot)
        self.SPiDist = 1. / self.SPDist        
        #make a shortcut array
        Y = np.zeros(self.NKnot+3, float)
        #fill in the y-values
        Y[1:-2] = self.Knots
        #make Y vals on end point cyclic
        Y[0] = self.Knots[-1]
        Y[-2] = self.Knots[0]
        Y[-1] = self.Knots[1]
        #make the coefficients
        C = np.zeros((4, self.NKnot), float)
        for i in range(1, self.NKnot + 1):
            C[:,i-1] = np.dot(self.BFunc, Y[i-1:i+3])
        self.SPC0 = C[0]
        self.SPC1 = C[1]
        self.SPC2 = C[2]
        self.SPC3 = C[3]
        
        
   
    
    