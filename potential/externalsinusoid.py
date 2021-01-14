#/usr/bin/env python
#Note: have not included derivatives wrt parameters yet

import base
import numpy as np


class ExternalSinusoid(base.PotentialClass):
    """External sinusoid field interactions."""

    Names = ["external_sinusoid"]   

    Source = """


#float PI
#PI = 4.d0*atan(1.0_8)
>>> defs
float PlaneDist
float ThisU
float factor

>>> argmainloopbeforepair
PlaneDist = (POSIMINIMAGE([[PlaneAxis]]) - PlaneLoc)
ARGVAL = PlaneDist
ARGTYPE = 0
!factor = 2*PI*NPeriods(0)/BOXL([[PlaneAxis]])
[ARGGET]

>>> argeval
!This section is for histogram evaluation
PlaneDist = ARGVAL
factor = 2*PI*NPeriods(0)/BOXL([[PlaneAxis]])
!THISU = UConst(0) * SIN(2*PI*(PlaneDist-PlaneLoc)*NPeriods(0)/BOXL([[PlaneAxis]]))
THISU = UConst(0) * SIN(factor*(PlaneDist-PlaneLoc))
!THISU = UConst(0) * PlaneDist
if (CalcDUParam) then
!    DU_UConst(0) = SIN(2*PI*(PlaneDist-PlaneLoc)*NPeriods(0)/BOXL([[PlaneAxis]]))
    DU_UConst(0) = SIN(factor*(PlaneDist-PlaneLoc))
endif
[UPDATEENERGY]

>>> mainloopbeforepair
!This section calculates the "explicit energy"
PlaneDist = (POSIMINIMAGE([[PlaneAxis]]) - PlaneLoc) 
if (BOXL([[PlaneAxis]])> 0.d0) PlaneDist = PlaneDist - BOXL([[PlaneAxis]]) * dnint(PlaneDist * IBOXL([[PlaneAxis]]))
factor = 2*PI*NPeriods(0)/BOXL([[PlaneAxis]])
!THISU = UConst(0) * SIN(2*PI*(PlaneDist-PlaneLoc)*NPeriods(0)/BOXL([[PlaneAxis]]))
THISU = UConst(0) * SIN(factor*(PlaneDist-PlaneLoc))
!THISU = UConst(0) * PlaneDist
[UPDATEENERGY]
if (CalcDUParam) then
!    DU_UConst(0) = SIN(2*PI*(PlaneDist-PlaneLoc)*NPeriods(0)/BOXL([[PlaneAxis]]))
    DU_UConst(0) = SIN(factor*(PlaneDist-PlaneLoc))
endif
if (CalcForce) then
    FORCEI = 0.d0
!    FORCEI([[PlaneAxis]]) = -UConst(0) * 2*PI*NPeriods(0)/BoxL([[PlaneAxis]]) * COS(2*PI*(PlaneDist-PlaneLoc)*NPeriods(0)/BOXL([[PlaneAxis]]))
    FORCEI([[PlaneAxis]]) = -UConst(0) * factor * COS(factor*(PlaneDist-PlaneLoc))
!    ForceI([[PlaneAxis]]) = - UConst(0)
Force(i,:) = Force(i,:) + Forcei
endif
""" 

    TestArgs = {"UConst":1.1, "PlaneAxis":1, "PlaneLoc":0.0}
    Type = base.ptypes.FieldPotential

    def __init__(self, Sys, Label = "", Filter = None,
                 Fixed = True,
                 UConst = None, NPeriods = 1.0, PlaneAxis = 0, PlaneLoc = 0.0):
        """Initializes an external sinusoid tied to current box size."""
        if PlaneAxis >= Sys.World.Dim:
            raise ValueError("PlaneAxis must be in range 0 to %d for this system." % (Sys.World.Dim-1))
        Source = ExternalSinusoid.Source.replace("[[PlaneAxis]]", "%d" % PlaneAxis)
        base.PotentialClass.__init__(self, Sys, Label = Label, 
                                     Filter = Filter, 
                                     Source = Source)
        self.ModuleVars += ["PlaneLoc", "PlaneAxis"]

        if UConst is None:
            raise ValueError("Must specify UConst.")
        
        self.Param.Add("UConst", 1, Value = UConst, Fixed = Fixed,
                       Scale = self.Sys.Units.EScale / self.Sys.Units.LScale)
        self.Param.Add("NPeriods", 1, Value = NPeriods, Fixed=Fixed,
                       Scale = self.Sys.Units.EScale / self.Sys.Units.LScale)
        self.PlaneLoc = PlaneLoc
        self.PlaneAxis = PlaneAxis
 
        self.Update()      
        

    def SetParam(self, UConst = None, NPeriods = None):
        """Sets parameters for this potential."""
        self.UConst = UConst
        self.NPeriods = NPeriods
        self.Update()

    def Val(self, x):
        """Returns the value of the potential."""
        return self.UConst[0] * np.sin( (x-self.PlaneLoc)*2*np.pi*self.NPeriods[0]/self.Sys.BoxL[self.PlaneAxis])
        
    def DVal(self, x ):
        """Returns the derivative of the potential."""
        return 2*np.pi*self.NPeriods[0]/self.Sys.BoxL[self.PlaneAxis] * self.UConst[0] * np.cos( (x-self.PlaneLoc)*2*np.pi*self.NPeriods[0]/self.Sys.BoxL[self.PlaneAxis])

    def Update(self):
        """Updates the potential."""
        pass

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.UConst.MaxChange = MaxFracChange * self.Sys.Units.EScale / self.Sys.Units.LScale

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        pass

    def LammpsStr(self):
        """Returns Appropriate Lammps String. """
        axisStrings=["lx","ly","lz"]
        
        #--- setting up variables ---
        s = "variable U{0} equal {1}\n".format(self.Label, self.UConst[0])
        s+= "variable offset{0} equal {1}\n".format(self.Label,self.PlaneLoc)
        s+= "variable nperiod{0} equal {1}\n".format(self.Label,self.NPeriods[0])
        s+= "variable factor{0} equal 2.0*PI*${{nperiod{0}}}/{1}\n".format(self.Label, axisStrings[self.PlaneAxis])
        s+= "variable sine{0} atom ${{U{0}}}*sin(${{factor{0}}}*(x-${{offset{0}}}))\n".format(self.Label)
        s+= "variable cosine{0} atom ${{U{0}}}*${{factor{0}}}*cos(${{factor{0}}}*(x-${{offset{0}}}))\n".format(self.Label)
        s+= "\n"

        #--- need to figure out types ---
        atomsInField = []
        #print(self.Filter)
        #print(self.Sys.World.SiteTypes)
        for Type in self.Filter.Select(self.Sys.World.SiteTypes):
            if type(Type) is list: #for some reason, if only one site type, Filter.Select returns a list instead of just that one item!
                Type = Type[0]
            
            #print(Type)
            #print(Type.SID)
            atomsInField.append("{:d}".format(Type.SID+1))
        
        typestr = " ".join(atomsInField)
        s+= "group atomsInField{0} type {1}\n".format(self.Label,typestr)
        
        #--- implement the force fix ---
        fxyz = ["0.0","0.0","0.0"]
        fxyz[self.PlaneAxis] = "v_cosine{0}".format(self.Label)
        s+= "fix Uext{0} atomsInField{0} addforce {1} energy v_sine{0}\n".format(self.Label," ".join(fxyz))
        s+= "fix_modify Uext{0} energy yes\n".format(self.Label)

        print(s)

        return s
                
            
        
        

        
