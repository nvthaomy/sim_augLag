#/usr/bin/env python


import numpy as np
import base


class Torsion(base.PotentialClass):
    """Torsion angle interactions, with cosine expansion."""

    Names = ["torsion"]    

    Source = """
>>> defs
int dind
float ang
float cosang
float sinang


>>> torsionenergydparam
ThisU = 0.d0
do dind = 0, [[NTerm - 1]]
    cosang =  cos(B(dind) * PHI + C(dind))
    ThisU = ThisU + A(dind) * (1.d0 + cosang)
enddo
[UPDATEENERGY]
if (CalcDUParam) then
    do dind = 0, [[NTerm - 1]]
        ang = B(dind) * PHI + C(dind)
        sinang = sin(ang)
        cosang = cos(ang)
        DU_A(dind) = 1.d0 + cosang
        DDU_A_B(dind,dind) = -PHI * sinang
        DDU_A_C(dind,dind) = -sinang
        DU_B(dind) = -A(dind) * PHI * sinang
        DDU_B_B(dind,dind) = -A(dind) * PHI*PHI * cosang
        DDU_B_C(dind,dind) = -A(dind) * PHI * cosang
        DU_C(dind) = -A(dind) * sinang
        DDU_C_C(dind,dind) = -A(dind) * cosang
    enddo
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
        ThisDU = 0.d0
        do dind = 0, [[NTerm - 1]]
            sinang = sin(B(dind) * PHI + C(dind))
            ThisDU = ThisDU - A(dind) * B(dind) * sinang
        enddo
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
    TestArgs = {"NTerm":2,
                "A":[1.2, 2.5],
                "B":[2.2, 3.7],
                "C":[3.1, 4.2]}
    Type = base.ptypes.TorsionPotential

    def __init__(self, Sys, Label = "", Filter = None,
                 Fixed = False,
                 NTerm = 1,
                 A = 0., B = 0., C = 0.):
        """Initializes a cosine torsion interaction."""
        Source = Torsion.Source.replace("[[NTerm - 1]]", "%d" % (NTerm - 1))
        base.PotentialClass.__init__(self, Sys, Label, Filter = Filter,
                                     Source = Source)

        self.NTerm = NTerm
        self.Param.Add("A", NTerm, Value = A, Fixed = Fixed,
                       Scale = self.Sys.Units.EScale)
        #test that B is integers
        l = B
        if not hasattr(B, "__getitem__"):
            l = [B]
        for x in l:
            if not abs(int(x) - x) < 1.e-10:
                print "Warning: B value %f is not an integer." % x
        self.Param.Add("B", NTerm, Value = B, Fixed = True, Scale = 1.)
        self.Param.Add("C", NTerm, Value = C, Fixed = Fixed, Scale = 1.)
#        self.Param.Add("C", NTerm, Value = C, Fixed = Fixed, Scale = 1., 
#                       Min = 0., SoftMin = 0.001,
#                       Max = 2.*np.pi, SoftMax = 1.99999 * np.pi)
        self.Update()


    def SetParam(self, A, B, C):
        """Sets parameters for this potential."""
        self.A = A
        self.B = B
        self.C = C
        self.Update()
        
    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.A.MaxChange = MaxFracChange * self.Sys.Units.EScale
            self.B.MaxChange = MaxFracChange * np.pi
            self.C.MaxChange = MaxFracChange * np.pi

    def Val(self, x):
        """Returns the value of the potential."""
        x = np.mod(x - self.C + np.pi, 2*np.pi) + self.C - np.pi
        return np.sum(self.A * (1. + np.cos(self.B*x + self.C)))
        
    def DVal(self, x):
        """Returns the value of the potential."""
        x = np.mod(x - self.C + np.pi, 2*np.pi) + self.C - np.pi
        return np.sum(-self.A * self.B * np.sin(self.B*x + self.C))        

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        self.SetParam(0, 0, 0)
        



        
        

        