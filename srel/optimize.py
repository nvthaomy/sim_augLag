#/usr/bin/env python


### Relative entropy optimizer -- general routines and class.
### coded by MSS


import numpy as np

import gaussjordan
import base
import penalty


#treatment for shiftable spline potentials
# 0 = nothing
# 1 = anchor knot with largest histogram to zero value
# 2 = constrain sum of knots to be zero 
SplineShiftOption = 1

#treatment for inner our outer regions in non-cyclic splines
# 0 = let edge knots vary
# 1 = let first two knots vary
SplineIndependentOption = 1

SplineVars = ["Knots", "Knots1", "Knots2"]
       
 
            
def GetBoundaryEnes(P):
    """Returns the target energies at the boundaries of a spline."""
    if not P.IsSpline:
        raise ValueError("Potential %s is not a spline potential" % P.Label)
    Units = P.Sys.Units
    Temp = P.Sys.TempSet
    EneInner = Units.ConvertEne(P.EneInner, Temp)
    EneSlopeInner = Units.ConvertEne(P.EneSlopeInner, Temp)
    EneOuter = Units.ConvertEne(P.EneOuter, Temp)
    EneSlopeOuter = Units.ConvertEne(P.EneSlopeOuter, Temp)
    EneInnerInit = Units.ConvertEne(P.EneInnerInit, Temp)
    EneOuterInit = Units.ConvertEne(P.EneOuterInit, Temp)
    EneSlopeInnerInit = Units.ConvertEne(P.EneSlopeInnerInit, Temp)
    EneSlopeOuterInit = Units.ConvertEne(P.EneSlopeOuterInit, Temp)
    EneSlopeInnerMin = Units.ConvertEne(P.EneSlopeInnerMin, Temp)
    EneSlopeOuterMin = Units.ConvertEne(P.EneSlopeOuterMin, Temp)    
    SlopeBuffer = Units.ConvertEne(P.EneSlopeBuffer, Temp)
    if P.KnotsCyclic and (EneInner != None or EneOuter != None or EneSlopeInner != None or EneSlopeOuter != None):
        raise ValueError("Cannot set potential slope or value at boundaries for cyclic potential %s" % P.Label)
    if P.KnotsShiftable and (EneInner != None or EneOuter != None):
        raise ValueError("Cannot set EneInner or EneOuter for potential %s because knots are shiftable.  Use EneSlopeInner / EneSlopeOuter instead." % P.Label)
    if EneOuter != None and EneSlopeOuter != None:
        raise ValueError("Cannot specify both the energy and slope at outer boundary for potential %s" % P.Label)
    if EneInner != None and EneSlopeInner != None:
        raise ValueError("Cannot specify both the energy and slope at inner boundary for potential %s" % P.Label)
    return (EneInner, EneSlopeInner, EneOuter, EneSlopeOuter, EneInnerInit, EneOuterInit, 
            EneSlopeInnerInit, EneSlopeOuterInit, EneSlopeInnerMin, EneSlopeOuterMin, SlopeBuffer)
    

class OptimizeClass(object):

    def __init__(self, ModSys,ElecSys):
        """Initializes a Srel-minimization class based on a reference target trajectory.
ModSys: System object for the model (coarse-grained) system; used for getting/setting constraints
"""
        if type(ModSys) is list:
            #filter for unique systems
            ModSys = [x for (i,x) in enumerate(ModSys) if not x in ModSys[i+1:]]
            #check for same parameters
            base.CheckAllModSys(ModSys)
            self.AllModSys = ModSys
            self.ModSys = ModSys[0]
        else:
            self.AllModSys = [ModSys]
            self.ModSys = ModSys
        if type(ElecSys) is list:
            #filter for unique systems
            ElecSys = [x for (i,x) in enumerate(ElecSys) if not x in ElecSys[i+1:]]
            self.AllElecSys = ElecSys
            self.ElecSys = ElecSys[0]
        else:
            self.AllElecSys = [ElecSys]
            self.ElecSys = ElecSys        
        self.Param = self.ModSys.ForceField.Param.Val
        self.Scale = self.ModSys.ForceField.Param.Scale
        self.Fixed = self.ModSys.ForceField.Param.Fixed
        self.Names = self.ModSys.ForceField.Param.Names
        self.dx = 0.
        #initialize constraints
        self.ClearConstraints()
        self.SplineKnotHists = None
        self.FLSplineKnotHists = None
        self.UseMaxChange = False
        #check for fluctuations
        self.HasFluct = False
        self.UseFluct = False
        for Sys in self.AllModSys:
            if any([P.UsesFluct for P in Sys.ForceField]):
                self.HasFluct = True
                self.UseFluct = True
                break
        self.Penalties = []
        self.UpdateMode = 0  #0 for normal, 1 for srel only, 2 for bias only
        self.Reset()        

    def __del__(self):
        """Garbage collecting routines."""
        pass
        

    def Reset(self):
        """Resets optimizer for completely new iterations."""
        self.Srel = 0.
        self.LastParam = None
        n = len(self.Param)
        self.LastSrel = None
        self.Bias = 0.
        self.FluctBias = 0.
        self.Obj = 0.
        self.DObj = np.zeros(n, float)
        self.DDObj = np.zeros((n,n), float)
        self.LastObj = None
        self.LastDObj = None
        self.Dir = np.zeros(n, float)
        self.LastDir = None
        self.SrelVals = []
        self.ObjVals = []
        self.BiasVals = []
        self.FluctBiasVals = []
        self.Iter = 0
        self.Mode = "INIT"
        
        
    def __ParseParamList(self, ParamList, Weights = None):
        """Gets parameter indices from a list of parameter objects."""
        NParam = len(self.Param)
        #make an array of weights
        if Weights is None:
            Weights = np.ones(len(ParamList), float)
        else:
            if len(Weights) == len(ParamList):
                Weights = np.array(Weights, dtype=float)
            else:
                raise TypeError("Weights must be same length as ParamList")
        #make a list of indices
        l1 = []
        w1 = []
        for (x,w) in zip(ParamList, Weights):
            if hasattr(x, "Ind1"):
                l1.extend(range(x.Ind1, x.Ind1 + len(x)))
                w1.extend([w]*len(x))
            elif type(x) in [int, np.dtype(int)]:
                if x < 0 or x >= NParam:
                    raise IndexError("%d in ParamList out of range." % x)
                else:
                    l1.append(int(x))
                    w1.append(w)
            else:
                raise TypeError("Unrecognized element %s in ParamList" % x)
        #get unique entries
        l2 = []
        w2 = []
        for (x,w) in zip(l1,w1):
            if x in l2:
                w2[l2.index(x)] += w
            else:
                l2.append(x)
                w2.append(w)
        #sort
        l3 = zip(l2, w2)
        l3.sort()
        l2 = [x for (x,w) in l3]
        w2 = [w for (x,w) in l3]
        return l2, w2
        
    def ClearConstraints(self):
        """Clears any parameter constraints."""
        self.Constraints = []
        self.ConstrParam = []
        self.Constrained = []
        self.ConstrIndep = []
        self.SplineConstraints = []
        self.Inequalities = []
                
    def ConstrainEquals(self, ParamList):
        """Constrains equality among a list of parameters.
ParamList: list of parameters (i.e., ModSys.ForceField.Param[2:4])"""
        NParam = len(self.Param)
        #get the indices of the supplied parameters
        l, w = self.__ParseParamList(ParamList)
        if len(l) < 2:
            raise ValueError("Parsed ParamList is too short.")
        #add constraints in the form of a constraint matrix
        for x in l[1:]:
            ConstrLine = np.zeros(NParam+1, float)
            ConstrLine[l[0]] = 1.
            ConstrLine[x] = -1.
            self.Constraints.append(ConstrLine)

    def ConstrainSum(self, ParamList, Value,
                     Weights = None):
        """Constrains sum of a list of parameters.
ParamList: list of parameters (i.e., ModSys.ForceField.Param[2:4])
Value: value of the sum of the parameters
Weights: optional list of weights"""
        NParam = len(self.Param)
        #get the indices of the supplied parameters
        l, Weights = self.__ParseParamList(ParamList, Weights)
        if len(l) < 1:
            raise ValueError("Parsed ParamList is too short.")
        #add constraints in the form of a constraint matrix
        ConstrLine = np.zeros(NParam+1, float)
        ConstrLine[-1] = Value
        for (i,x) in enumerate(l):
            ConstrLine[x] = Weights[i]
        self.Constraints.append(ConstrLine)
        
    def ConstrainGreaterThan(self, ParamList, Min, SoftMin,
                             Weights = None):
        """Constrains sum of a list of parameters to be greater than a number.
ParamList: list of parameters (i.e., ModSys.ForceField.Param[2:4])
Min: value of the sum of the parameters to be greater than
SoftMin: when the penalty turns on (then diverges at Value)
Weights: optional list of weights"""
        NParam = len(self.Param)
        #get the indices of the supplied parameters
        l, Weights = self.__ParseParamList(ParamList, Weights)
        if len(l) < 1:
            raise ValueError("Parsed ParamList is too short.")
        #add constraints in the form of a constraint matrix
        ConstrLine = np.zeros(NParam, float)
        for (i,x) in enumerate(l):
            ConstrLine[x] = Weights[i]
        if SoftMin <= Min:
            raise ValueError("SoftMin must be greater than Min in inequality constraint.")
        self.Inequalities.append((ConstrLine, SoftMin, Min))    
        
    def ConstrainChargePerAtom(self, ChargePerAtom = 0.):
        """Constrains the overall charge per atom.
ChargePerAtom: total charge per atom for the entire system (dflt 0.)"""
        Charge = self.ModSys.ForceField.Globals.Charge
        NAtom = self.ModSys.AIDCount.astype(float)
        ParamList = [Charge.Ind1 + i for i in range(len(Charge))]
        ChargeSum = np.dot(Charge, NAtom)
        self.ConstrainSum(ParamList, ChargeSum)
        
    def ConstrainNeutralCharge(self, ChargePerAtom = 0.):
        """Ensures the overall system charge is zero."""
        Charge = self.ModSys.ForceField.Globals.Charge
        Weights = self.ModSys.AIDCount.astype(float)
        ParamList = [Charge.Ind1 + i for i in range(len(Charge))]
        self.ConstrainSum(ParamList, 0., Weights)        


    def ConstrainSplines(self): 
        """Constrains any spline knots that are not sampled based on various options.
"""
        NConstraints = len(self.Constraints)
        
        #loop through potentials 
        for (ind, P) in enumerate(self.ModSys.ForceField):
            
            #is this a spline?
            if not P.IsSpline: continue
            
            #should we constrain this spline?
            if not P.ConstrainSpline: continue
           
            #get the boundary options for this spline
            (EneInner, EneSlopeInner, EneOuter, EneSlopeOuter, EneInnerInit, EneOuterInit,
             EneSlopeInnerInit, EneSlopeOuterInit, EneSlopeInnerMin, EneSlopeOuterMin, SlopeBuffer) = GetBoundaryEnes(P)
            DistPerKnot = P.SPDist

            #get the histogram for this potential
            if self.SplineKnotHists is None:
                KnotHist = P.KnotHist()
            else:
                KnotHist = self.SplineKnotHists[ind]    
                
            #set the spline options
            Name = P.Name
            self.__ConstrainKnots(Name, P.Knots, KnotHist, DistPerKnot,
                P.ForceInner, P.ForceOuter, P.KnotsCyclic, P.KnotsShiftable,
                P.KnotMinHistFrac, P.KnotMinHistFracInner, P.KnotMinHistFracOuter,
                ValInner = EneInner, SlopeInner = EneSlopeInner, 
                InitInner = EneInnerInit, InitSlopeInner = EneSlopeInnerInit,
                ValOuter = EneOuter, SlopeOuter = EneSlopeOuter, 
                InitOuter = P.EneInnerInit, InitSlopeOuter = EneSlopeOuterInit,
                MinSlopeInner = EneSlopeInnerMin, MinSlopeOuter = EneSlopeOuterMin,
                SlopeBuffer = SlopeBuffer)
            
            #now constrain the energy fluctuation knots
            if P.UsesFluct:
                #get the histogram for this potential
                if self.FLSplineKnotHists is None:
                    KnotHist = P.FLKnotHist()
                else:
                    KnotHist = self.FLSplineKnotHists[ind]  
                #set the spline options
                DistPerKnot = P.FLDist
                if P.FluctTreatment == 0:
                    #make fluctuations zero in unsampled regions
                    ValInner = P.FluctSoftMin
                    ValOuter = P.FluctSoftMin
                    InitInner = P.FluctSoftMin
                    InitOuter = P.FluctSoftMin
                    self.__ConstrainKnots(Name, P.FLKnots, KnotHist, DistPerKnot,
                        False, False, P.KnotsCyclic, False,
                        P.KnotMinHistFrac, P.KnotMinHistFracInner, P.KnotMinHistFracOuter,
                        ValInner=ValInner, InitInner=InitInner,
                        ValOuter=ValOuter, InitOuter=InitOuter)
                else:
                    #make slope of fluctuations zero in unsampled regions
                    self.__ConstrainKnots(Name, P.FLKnots, KnotHist, DistPerKnot,
                        False, False, P.KnotsCyclic, False,
                        P.KnotMinHistFrac, P.KnotMinHistFracInner, P.KnotMinHistFracOuter,
                        SlopeInner = 0., SlopeOuter = 0.)
                                     
        #separate out the spline constraints
        self.SplineConstraints = self.Constraints[NConstraints:]
        self.Constraints = self.Constraints[:NConstraints]
            

    def __ConstrainKnots(self, Name, Knots, KnotHist, DistPerKnot,
                         ForceInner, ForceOuter, KnotsCyclic, KnotsShiftable,
                         MinHistFrac, MinHistFracInner, MinHistFracOuter,
                         ValInner=None, SlopeInner=None, InitInner=None, InitSlopeInner=None,
                         ValOuter=None, SlopeOuter=None, InitOuter=None, InitSlopeOuter=None,
                         MinSlopeInner = None, MinSlopeOuter = None, SlopeBuffer = None):
        """Helper routine for doing spline constraints for parameters in Knots."""   
            
        #number of knots
        n = len(Knots)
        #starting index in parameter array of knots
        Shift = Knots.Ind1
    
        #find the minimum value of the histogram below which corresp. knots will be constrained
        MinHist = MinHistFrac * np.max(KnotHist)
        MinHistInner = MinHistFracInner * np.max(KnotHist)
        MinHistOuter = MinHistFracOuter * np.max(KnotHist)
        #check in case the maximum is zero, for which this potential isn't used at all so skip 
        if np.max(KnotHist) == 0: 
            return
        if base.DEBUG: 
            print "Spline histograms for potential %s" % Name
            print "  Minimum histogram value: %.3e" % MinHist
            print "  Minimum inner edge histogram value: %.3e" % MinHistInner
            print "  Minimum outer edge histogram value: %.3e" % MinHistOuter
            for (i, h) in enumerate(KnotHist):
                print "  %s histogram value: %.3e" % (self.Names[i+Shift], h)
                
        #make groups of knots that are unsampled
        KnotGroups = []
        KnotsUnsampled = set()
        
        #special treatment for the inner and outer regions
        if not KnotsCyclic:
            KnotsI = set()
            KnotsO = set()
            #examine inner region; add all knots until there are samples
            for (i, h) in enumerate(KnotHist):
                if h >= MinHistInner:
                    break
                KnotsUnsampled.add(i)
            if i > 0 and i < n - 2:
                KnotsI = set(range(0, i+1))
            #examine outer region; add all knots up to where there are samples
            for (i, h) in reversed(list(enumerate(KnotHist))):
                if h >= MinHistOuter:
                    break
                KnotsUnsampled.add(i)
            if i > 2 and i < n-1:
                KnotsO = set(range(i,n))
            #check for overlap
            if KnotsI & KnotsO:
                midpoint = np.mean(list(KnotsI & KnotsO))
                KnotsI = set([x for x in KnotsI if x <= midpoint])
                KnotsO = set([x for x in KnotsO if x > midpoint])
            if KnotsI: KnotGroups.append(KnotsI)
            if KnotsO: KnotGroups.append(KnotsO)
        
        #examine other regions, not inner or outer
        for (i, h) in enumerate(KnotHist):
            if h <= MinHist:
                KnotsUnsampled.add(i)
                if i == 0:
                    if KnotsCyclic:
                        Knots = set([n-1, 0, 1])
                    else:
                        Knots = set([0, 1])
                elif i == n-1:
                    if KnotsCyclic:
                        Knots = set([n-2, n-1, 0])
                    else:
                        Knots = set([n-2, n-1])
                else:
                    Knots = set([i-1, i, i+1])
                KnotGroups.append(Knots)

        #add knots according to boundary conditions
        if not KnotsCyclic:
            DoInner = (0 in KnotsUnsampled) or ForceInner
            if DoInner and not (ValInner is None and SlopeInner is None):
                KnotGroups.append(set([0, 1]))
            DoOuter = (n-1 in KnotsUnsampled) or ForceOuter
            if DoOuter and not (ValOuter is None and SlopeOuter is None):
                KnotGroups.append(set([n-2, n-1]))
            
        #now make super sets out of these by combining groups with overlap
        #each group will involve neighbors on either end which are not constr.
        Cont = True
        while Cont:
            Cont = False
            i = 0
            j = 1
            while i < len(KnotGroups) and j < len(KnotGroups):
                Knotsi = KnotGroups[i]
                Knotsj = KnotGroups[j]
                if len(Knotsi & Knotsj) > 0:
                    Knotsij = Knotsi | Knotsj
                    KnotGroups[i] = Knotsij
                    del KnotGroups[j]
                    Cont = True
                else:
                    j += 1
                    if j >= len(KnotGroups):
                        i += 1
                        j = i+1
                        
        #make sorted lists of all unique supersets of constraint groups and add base index;
        #special care to treat cyclic knots
        for (indKnots, Knots) in enumerate(KnotGroups): 
            Knots = sorted(list(Knots))
            if 0 in Knots and n-1 in Knots and KnotsCyclic:
                #must do this in a special way to ensure correct first and last knots
                #for example, want to transform knots [0, 1, 2, 18, 19] to [18, 19, 0, 1, 2]
                for i in Knots:
                    if not i+1 in Knots and not i == n-1: iMax = i
                    if not i-1 in Knots and not i == 0: iMin = i
                #remake the knot group
                Knots = [i for i in Knots if i >= iMin] + [i for i in Knots if i <= iMax]
                if base.DEBUG:
                    print "Found cyclic group of knots and reordered:"
                    print Knots
            #add base index
            Knots = [Shift + x for x in Knots]
            #update knot group
            KnotGroups[indKnots] = Knots
            
        FirstKnot = Shift
        LastKnot = Shift + n - 1            
            
        #can we shift this potential?
        CanShift = KnotsShiftable and (ValInner == None) and (ValOuter == None)

        #now go through and add constraints
        for Knots in KnotGroups:
            
            Thisn = len(Knots)
            
            #check if need to constrain this group of knots
            if FirstKnot in Knots and not KnotsCyclic:
                #constrain parameters in inner region
                iMax = max(Knots)
                if iMax == LastKnot:
                    raise ValueError("First and last knots in potential %s are constrained." % Name)
                if not SlopeInner is None:
                    #constrain parameters to be linear with specified slope
                    DKnot = SlopeInner * DistPerKnot
                    for i in reversed(Knots[1:]):
                        self.ConstrainSum([i-1, i], DKnot, Weights = [1, -1])
                        self.Param[i-1] = self.Param[i] + DKnot
                    #prioritize last knot in group as independent
                    self.ConstrIndep.extend([iMax])
                    if base.DEBUG:
                        print "Spline inner linear ramp constraint group, slope=%f:" % SlopeInner                        
                elif not ValInner is None:
                    self.Param[FirstKnot] = ValInner
                    #constrain value at inner edge with linear slope up until then
                    DKnot = -(ValInner - self.Param[iMax]) / (Thisn - 1.)
                    self.ConstrainSum([FirstKnot], ValInner)
                    for i in Knots[1:-1]:
                        #constrain Knots[i] - Knots[i-1] equal to (Knots[iMax] - Knots[FirstKnot]) / (iMax - FirstKnot)
                        f = 1. / float(iMax - FirstKnot)
                        self.ConstrainSum([i, i-1, iMax, FirstKnot], 0., 
                                          Weights = [-1, 1, f, -f])
                        self.Param[i] = self.Param[i-1] + DKnot
                    #prioritize last knot in group as independent
                    self.ConstrIndep.extend([iMax])
                    if base.DEBUG:
                        print "Spline inner linear ramp constraint group, value=%f:" % ValInner
                else:
                    #constrain parameters to be an extrapolated linear ramp
                    if not InitInner is None:
                        self.Param[FirstKnot] = InitInner
                    elif not InitSlopeInner is None:
                        self.Param[FirstKnot] = self.Param[iMax] + InitSlopeInner * DistPerKnot * (Thisn - 1.)
                    DKnot = -(self.Param[FirstKnot] - self.Param[iMax]) / (Thisn - 1.)
                    for i in Knots[1:-1]:
                        #constraint Knots[i] - Knots[i-1] equal to Knots[i1] - Knots[i1-1]
                        self.ConstrainSum([i, i-1, iMax, iMax-1], 0., 
                                          Weights = [-1, 1, 1, -1])
                        self.Param[i] = self.Param[i-1] + DKnot
                    if SplineIndependentOption == 0:
                        #prioritize first and last knots in group as independent
                        self.ConstrIndep.extend([FirstKnot, iMax])
                    else:
                        #prioritize penultimate and last knots in group as independent
                        self.ConstrIndep.extend([iMax-1, iMax])
                    if not MinSlopeInner is None:
                        #add an inequality constraint so the first knot is at least buffer more than the second knot
                        MinDelta = DistPerKnot * MinSlopeInner * (Thisn - 1.)
                        self.ConstrainGreaterThan([Knots[0], Knots[-1]], Min = MinDelta - SlopeBuffer,
                                                  SoftMin = MinDelta, Weights = [1., -1])  
                    if base.DEBUG:
                        print "Spline inner linear ramp constraint group:"
                if base.DEBUG:
                    for i in Knots:
                        print "  %s" % self.Names[i]
                        
            elif LastKnot in Knots and not KnotsCyclic:
                #constrain parameters in outer region
                iMin = min(Knots)
                if iMin == FirstKnot:
                    raise ValueError("First and last knots in potential %s are constrained." % Name)
                if not SlopeOuter is None:
                    #constrain parameters to be linear with specified slope
                    DKnot = SlopeOuter * DistPerKnot
                    for i in Knots[1:]:
                        self.ConstrainSum([i-1, i], DKnot, Weights = [-1, 1]) 
                        self.Param[i] = self.Param[i-1] + DKnot
                    #prioritize first knot in group as independent
                    self.ConstrIndep.extend([iMin])  
                    if base.DEBUG:
                        print "Spline outer linear ramp constraint group, slope=%f:" % SlopeOuter
                elif not ValOuter is None:
                    #constrain value edge, with linear slope up until then
                    self.Param[LastKnot] = ValOuter
                    DKnot = (ValOuter - self.Param[iMin]) / (Thisn - 1.)
                    self.ConstrainSum([LastKnot], ValOuter)
                    for i in Knots[1:-1]:
                        #constraint Knots[i] - Knots[i-1] equal to (Knots[LastKnot] - Knots[iMin]) / (LastKnot - iMin)
                        f = 1. / float(LastKnot - iMin)
                        self.ConstrainSum([i, i-1, LastKnot, iMin], 0., 
                                          Weights = [-1, 1, f, -f])  
                        self.Param[i] = self.Param[i-1] + DKnot
                    #prioritize last knot in group as independent
                    self.ConstrIndep.extend([iMin])    
                    if base.DEBUG:
                        print "Spline outer linear ramp constraint group, value=%f:" % ValOuter
                else:
                    #constrain parameters to be an extrapolated linear ramp
                    if not InitOuter is None:
                        self.Param[LastKnot] = InitOuter 
                    elif not InitSlopeOuter is None:
                        self.Param[LastKnot] = self.Param[iMin] + InitSlopeOuter * DistPerKnot * (Thisn - 1.)                        
                    DKnot = (self.Param[LastKnot] - self.Param[iMin]) / (Thisn - 1.)
                    for i in Knots[1:-1]:
                        #constraint Knots[i+1] - Knots[i] equal to Knots[iMin+1] - Knots[iMin]
                        self.ConstrainSum([i+1, i, iMin+1, iMin], 0., 
                                          Weights = [-1, 1, 1, -1])
                        self.Param[i] = self.Param[i-1] + DKnot
                    if SplineIndependentOption == 0:
                        #prioritize first and last knot in group as independent
                        self.ConstrIndep.extend([LastKnot, iMin])
                    else:
                        #prioritize first and second knots as independent
                        self.ConstrIndep.extend([iMin+1, iMin])
                    if not MinSlopeOuter is None:
                        #add an inequality constraint so the first knot is at least buffer more than the second knot
                        MinDelta = DistPerKnot * MinSlopeOuter * (Thisn - 1.)
                        self.ConstrainGreaterThan([Knots[0], Knots[-1]], Min = MinDelta - SlopeBuffer, 
                                                  SoftMin = MinDelta, Weights = [-1., 1.])                         
                    if base.DEBUG:
                        print "Spline outer linear ramp constraint group:"
                if base.DEBUG:
                    for i in Knots:
                        print "  %s" % self.Names[i]
                    
            else:
                #constrain parameters to be a linear ramp
                #this assumes that knots are properly ordered
                iFirst = Knots[0]
                iFirstp1 = Knots[1]
                iLast = Knots[-1]
                DKnot = (self.Param[iLast] - self.Param[iFirst]) / (Thisn - 1.)
                indLast = len(Knots) - 1
                for (indKnot, i) in enumerate(Knots):
                    #skip first and last knots
                    if indKnot == 0 or indKnot == indLast: continue
                    #find the previous and next knots
                    im1 = Knots[indKnot - 1]
                    ip1 = Knots[indKnot + 1]
                    #constrain Knots[indKnot+1] - Knots[indKnot] equal to Knots[1] - Knots[0]
                    self.ConstrainSum([ip1, i, iFirstp1, iFirst], 0., 
                                      Weights = [-1, 1, 1, -1])
                    self.Param[i] = self.Param[im1] + DKnot
                if base.DEBUG:
                    print "Spline linear ramp constraint group:"
                    for i in Knots:
                        print "  %s" % self.Names[i]
                #prioritize first and last knost in group as independent
                self.ConstrIndep.extend([iFirst, iLast])
                                                    
        #constrain sum of knots to be zero in cyclic splines
        if KnotsShiftable and CanShift:
            if SplineShiftOption == 1:
                i = np.argmax(KnotHist)
                self.ConstrainSum([i + Shift], 0.)
            elif SplineShiftOption == 2:
                self.ConstrainSum(range(Shift, Shift + n), 0.) 
                Knots = Knots - np.mean(Knots)
                                 

    def InitConstraints(self, Verbose = False):
        """Prepares the list of linear constraints."""
        self.ConstrainSplines()
        Constraints = np.array(self.Constraints + self.SplineConstraints, float)
        if len(Constraints):
            ConstrList = gaussjordan.DependentVars(Constraints, IndepPriorities = self.ConstrIndep)
        else:
            ConstrList = []
        self.ConstrParam = []
        self.Constrained = []
        #gow through the constraints
        for (Ind, Coef, Const) in ConstrList:
            #only add a constraint if not already fixed
            self.ConstrParam.append((Ind, Coef, Const))
            self.Constrained.append(Ind) 
            if base.DEBUG:
                print "Constrained param %s" % self.Names[Ind]
                print "      const: %f" % Const
                for (n,c) in zip(self.Names, Coef):
                    if abs(c) < 1.e-200: continue
                    print "  %s coef: %f" % (n, c)  
        for (Coef, SoftMin, Min) in self.Inequalities:
            #only add a constraint if not already fixed
            if base.DEBUG:
                print "Inequality constraint:"
                print "   soft min: %f" % Min
                print "        min: %f" % Min
                for (n,c) in zip(self.Names, Coef):
                    if abs(c) < 1.e-200: continue
                    print "  %s coef: %f" % (n, c)                         
        #initialize values
        self.SetParam(self.Param)
        if Verbose and len(self.ConstrParam) > 0:
            print "\nLinearly constrained parameters:"
            for (Ind, Coef, Const) in self.ConstrParam:
                if not self.Fixed[Ind]:
                    print "  " + self.Names[Ind]
            print ""
            
            
    
    def AddPenalty(self, Measure, Target, Coef = 1., Name = None, ValInd = 0, MeasureScale = 1.):
        """Adds a penalty constraint in the Srel minimization, of the form
Coef * (<Measure> - Target)^2  where <Target> is the ensemble average from Measure.  
Coef can be gradually increased over successive minimizations until the constraint
is satisfied within reasonable tolerance.
Measure:  An instance of a MeasureAvgClass object that is used to compute 
          the average to be constrained, or strings "PEnergy" or "Pressure"
Target:  The target value of the average
Coef:  The coefficient of the penalty function when added to the main objective
Name:  An optional name
ValInd: Index of the value to use in measure, if more than one
MeasureScale: optional scale factor for the measurement and target value""" 
        if Name is None:
            Name = "penalty%d" % len(self.Penalties)
        Pen = penalty.PenaltyClass(self.ModSys, Measure, Target, Coef, 
                                   Name = Name, ValInd = ValInd, MeasureScale = MeasureScale)     
        self.Penalties.append(Pen)
        return Pen
        
    def DelPenalty(self, Pen):
        """Removes a penalty constraint from the minimization."""
        if not Pen in self.Penalties:
            raise ValueError("Cound not find penalty in current list.")
        self.Penalties = [x for x in self.Penalties if not x is Pen]
        
                   
    def SetParam(self, Param):
        """Sets parameters, using constraints."""
        base.CheckAllModSys(self.AllModSys)
        #make a copy of parameters
        Param = Param.copy()
        #apply constraints
        for (Ind, Coef, Const) in self.ConstrParam:
            Param[Ind] = np.dot(Param, Coef) + Const
        #check min and max vals
        Min = self.ModSys.ForceField.Param.Min
        Max = self.ModSys.ForceField.Param.Max
        Violated = np.nonzero(np.logical_or(Param > Max, Param < Min))[0]
        if len(Violated) > 0:
            s = "The following parameters were set out of Min, Max range:\n"
            for i in Violated:
                s += "   parameter %d (%s)\n" % (i, self.Names[i])
            raise ValueError(s)
        #update systems with current parameters
        for i,ModSys in enumerate(self.AllModSys):
            ModSys.ForceField.SetParam(Param, ForceFixed = True)
            s = self.ModSys.ForceField.ParamString()
#            print('ModSys Param:\n{}'.format(s))
            ElecSys = self.AllElecSys[i]
            if ElecSys != None:
                ElecSys.ForceField.SetParamString(s)
#                print('ElecSys Param:\n{}'.format(ElecSys.ForceField.ParamString()))

    def CalcDistLim(self):
        """Computes the distance along direction Dir until a bound is reached
or until a maximum change in a parameter is reached.
Returns: distance until the hard limit, distance until the soft limit"""
        #find which variable will exceed min, max range first
        FF = self.ModSys.ForceField
        Min = FF.Param.Min
        Max = FF.Param.Max
        DistToMin = np.zeros_like(self.Param)
        DistToMax = np.zeros_like(self.Param)
        #update direction based on constraints
        Dir = self.Dir.copy()
        for (Ind, Coef, Const) in self.ConstrParam:
            Dir[Ind] = np.dot(Dir, Coef)
        #find distances to min and max
        for (i, d) in enumerate(Dir):
            if d == 0: continue
            try:
                DistToMin[i] = (Min[i] - self.Param[i]) / d
            except FloatingPointError:
                DistToMin[i] = 1.e300
            try:
                DistToMax[i] = (Max[i] - self.Param[i]) / d
            except FloatingPointError:
                DistToMax[i] = 1.e300
        DistToMin[DistToMin <= 0] = 1.e300
        DistToMax[DistToMax <= 0] = 1.e300
        DistTo = np.where(DistToMin < DistToMax, DistToMin, DistToMax)
        DistLim = np.min(DistTo)
        #find distances to soft min and max
        SoftMin = FF.Param.SoftMin
        SoftMax = FF.Param.SoftMax
        DistToMin = np.zeros_like(self.Param)
        DistToMax = np.zeros_like(self.Param)
        for (i, d) in enumerate(Dir):
            if d == 0: continue
            try:
                DistToMin[i] = (SoftMin[i] - self.Param[i]) / d
            except FloatingPointError:
                DistToMin[i] = 1.e300
            try:
                DistToMax[i] = (SoftMax[i] - self.Param[i]) / d
            except FloatingPointError:
                DistToMax[i] = 1.e300
        DistToMin[DistToMin <= 0] = 1.e300
        DistToMax[DistToMax <= 0] = 1.e300
        DistTo = np.where(DistToMin < DistToMax, DistToMin, DistToMax)
        SoftDistLim = np.min(DistTo)
        #add inequality constraints minimum
        for (i, (Coef, SoftMin, Min)) in enumerate(self.Inequalities):
            try:
                DistToMin = (Min - np.dot(self.Param, Coef)) / np.dot(Dir, Coef)
            except FloatingPointError:
                DistToMin = 1.e300
            if DistToMin > 0: 
                DistLim = min(DistLim, DistToMin)
        #add inequality constraints soft minimum
        for (i, (Coef, SoftMin, Min)) in enumerate(self.Inequalities):
            try:
                DistToMin = (SoftMin - np.dot(self.Param, Coef)) / np.dot(Dir, Coef)
            except FloatingPointError:
                DistToMin = 1.e300
            if DistToMin > 0: 
                SoftDistLim = min(SoftDistLim, DistToMin) 
        #now find which variable will exceed maxchange limit first
        if self.UseMaxChange:
            MaxChange = FF.Param.MaxChange
            DistToMaxChange = np.zeros_like(self.Param)
            DistToMaxChange[self.Dir!=0] = np.abs(MaxChange[self.Dir!=0] / self.Dir[self.Dir!=0])
            DistToMaxChange[self.Dir==0] = 1.e300
            DistLimMaxChange = np.min(DistToMaxChange)
        else:
            DistLimMaxChange = 1.e300
        #now combine range limits and max change limits
        if DistLimMaxChange < SoftDistLim:
            return DistLimMaxChange, 0.5 * DistLimMaxChange
        elif DistLimMaxChange < DistLim:
            return DistLimMaxChange, SoftDistLim
        else:
            return DistLim, SoftDistLim

    def PrepUpdate(self):
        """Records last parameters in preparation for an iteration."""
        self.LastObj = self.Obj
        self.LastDObj = self.DObj.copy()
        self.LastSrel = self.Srel
        self.LastParam = self.Param.copy()


    def UpdateObj(self):
        """Updates the objective function, adjusting for biases, penalties, and constraints.
This should be called after Obj, DObj, and DDObj are already filled with pure
relative entropy contributsions."""
        #The bias designed to keep all parameters within specified range from
        #range Min to Max is as follows: it is zero between SoftMin and SoftMax
        #and increases sharply, to infinity, from SoftMin to Min and
        #from SoftMax to Max.  The form of the increase is u^3 / (1-u)
        #where u is (x-SoftMin)/(Min-SoftMin) or (x-SoftMax)/(Max-SoftMax)
        #and where x is the value of the parameter.  The bias is also
        #scaled by the number of atoms in the system times 1000 to be
        #roughly on the same order of size as Srel.
        #
        #UpdateMode = 0 for full objective, 1 for srel only, 2 for bias only

        #initialize things
        if self.UpdateMode == 2:
            self.Obj = 0.
            self.DObj.fill(0.)
            self.DDObj.fill(0.)
        
        #shortcuts
        FF = self.ModSys.ForceField
        Min = FF.Param.Min
        SoftMin = FF.Param.SoftMin
        Max = FF.Param.Max
        SoftMax = FF.Param.SoftMax
        
        if not self.UpdateMode == 1:
            
            #scale range bias by system size
            BiasScale = float(len(self.ModSys.Atom)) * 1000
            #loop through parameters and add range bias
            for i in range(FF.Param.N):
                if self.Fixed[i]: continue
                if self.Param[i] < SoftMin[i]:
                    a = SoftMin[i]
                    b = Min[i]
                    l = 1./(b - a)
                    u = (self.Param[i] - a) * l
                    ThisBias = BiasScale * u*u*u / (1 - u)
                    self.Bias += ThisBias
                    self.Obj += ThisBias
                    self.DObj[i] += BiasScale * (u*u*(3-2*u)/(1-u)**2) * l
                    self.DDObj[i,i] += BiasScale * (2*u*(3-3*u+u*u)/(1-u)**3) * l*l
                elif self.Param[i] > SoftMax[i]:
                    a = SoftMax[i]
                    b = Max[i]
                    l = 1./(b - a)
                    u = (self.Param[i] - a) * l
                    ThisBias = BiasScale * u*u*u / (1 - u)
                    self.Bias += ThisBias
                    self.Obj += ThisBias
                    self.DObj[i] += BiasScale * (u*u*(3-2*u)/(1-u)**2) * l
                    self.DDObj[i,i] += BiasScale * (2*u*(3-3*u+u*u)/(1-u)**3) * l*l
                    
            #loop through inequality constraints
            for (Coef, ThisSoftMin, ThisMin) in self.Inequalities:
                Val = np.dot(Coef, self.Param)                
                if Val < ThisSoftMin:
                    a = ThisSoftMin
                    b = ThisMin
                    l = 1./(b - a)
                    u = (Val - a) * l
                    ThisBias = BiasScale * u*u*u / (1 - u)
                    self.Bias += ThisBias
                    self.Obj += ThisBias
                    self.DObj += BiasScale * Coef * (u*u*(3-2*u)/(1-u)**2) * l
                    self.DDObj += BiasScale * np.outer(Coef, Coef) * (2*u*(3-3*u+u*u)/(1-u)**3) * l*l                  
    
            #now add the penalties
            for Pen in self.Penalties:
                self.Obj, self.Bias, self.DObj, self.DDObj = Pen.UpdateObj(self.Obj, self.Bias, 
                                                                           self.DObj, self.DDObj)
    
        #update lists
        self.SrelVals.append(self.Srel)
        self.ObjVals.append(self.Obj)
        self.BiasVals.append(self.Bias)
        self.FluctBiasVals.append(self.FluctBias)
        self.Iter = len(self.ObjVals)
        
        #update derivative based on constraints
        NParam = len(self.Param)
        #make a matrix for linear mapping of constrained params
        A = np.identity(NParam, dtype=float)
        for (Ind1, Coef1, Const1) in self.ConstrParam:
            A[Ind1,:] = Coef1
        #update the derivatives
        self.DObj = np.dot(self.DObj.T, A) 
        self.DDObj = np.dot(A.T, np.dot(self.DDObj, A))


    def Output0(self):
        s = "INITIAL PARAMETERS\n"
        s += "   %-20s %-12s\n" % ("Label", "Param")
        for (i, nm) in enumerate(self.Names):
            if self.Fixed[i]: continue
            s += "   %-20s %-12.4e\n" % (nm, self.Param[i])
        return s 

    def Output1(self, OtherData = [], ShowParams = True):
        if hasattr(self, "ReweightFrac"):
            OtherData = [("ReweightFrac", self.ReweightFrac)] + OtherData
        s = "ITERATION %d  (%s)\n" % (self.Iter, self.Mode)
        if len(self.SrelVals) < 2:
            FracSrel = 0.
            FracObj = 0.
            FracBias = 0.
            FracFluctBias = 0.
        else:
            FracSrel = (self.SrelVals[-1] - self.SrelVals[-2]) / (abs(self.SrelVals[-1]) + 1.e-300)
            FracObj = (self.ObjVals[-1] - self.ObjVals[-2]) / (abs(self.ObjVals[-1]) + 1.e-300)
            FracBias = (self.BiasVals[-1] - self.BiasVals[-2]) / (abs(self.BiasVals[-1]) + 1.e-300)
            FracFluctBias = (self.FluctBiasVals[-1] - self.FluctBiasVals[-2]) / (abs(self.FluctBiasVals[-1]) + 1.e-300)
        Frac = (self.Param - self.LastParam) / (abs(self.LastParam) + 1.e-300) * 100.
        s += "   %-20s %12s %9s\n"     % ("Variable", "Value", "FracChnge")
        s += "   %-20s %12.4e %9.1e\n" % ("Obj", self.Obj, FracObj)
        s += "   %-20s %12.4e %9.1e\n" % ("Srel", self.Srel, FracSrel)
        s += "   %-20s %12.4e %9.1e\n" % ("Bias", self.Bias, FracBias)
        if self.HasFluct:
            s += "   %-20s %12.4e %9.1e\n" % ("FluctBias", self.FluctBias, FracFluctBias)        
        s += "   %-20s %12.4e\n"       % ("dx", self.dx)
        for (k,v) in OtherData:
            if type(v) in [np.ndarray, list]:
                for (i,x) in enumerate(v):
                    s += "   %-20s %12.4e\n" % (k.strip() + str(i), x)
            else:
                s += "   %-20s %12.4e\n" % (k, v)
        for Pen in self.Penalties:
            s += "   %-20s %12.4e\n" % (Pen.Name + "_avg", Pen.Avg)
        if not ShowParams:
            return s
        for (i, nm) in enumerate(self.Names):
            if self.Fixed[i]: continue
            s += "   %-20s %12.4e %9.1e\n" % (nm, self.Param[i], Frac[i])
        return s

    def Output2(self, OtherData = [], ShowParams = True):
        if hasattr(self, "ReweightFrac"):
            OtherData = [("ReweightFrac", self.ReweightFrac)] + OtherData
        s = "ITERATION %d  (%s)\n" % (self.Iter, self.Mode)
        s += "   %-20s %12s %12s %12s %12s\n" % ("Variable", "Value",
                                                 "DObj", "DDObj", "SearchDir")
        s += "   %-20s %12.4e\n" % ("Obj", self.Obj)
        s += "   %-20s %12.4e\n" % ("Srel", self.Srel)
        s += "   %-20s %12.4e\n" % ("Bias", self.Bias)
        if self.HasFluct:
            s += "   %-20s %12.4e\n" % ("FluctBias", self.FluctBias)   
        s += "   %-20s %12.4e\n" % ("dx", self.dx)
        for (k,v) in OtherData:
            if type(v) in [np.ndarray, list]:
                for (i,x) in enumerate(v):
                    s += "   %-20s %12.4e\n" % (k.strip() + str(i), x)
            else:
                s += "   %-20s %12.4e\n" % (k, v)
        for Pen in self.Penalties:
            s += "   %-20s %12.4e\n" % (Pen.Name + "_avg", Pen.Avg)
        if not ShowParams:
            return s                
        for (i, nm) in enumerate(self.Names):
            if self.Fixed[i]: continue
            s += "   %-20s %12.4e %12.4e %12.4e %12.4e\n" % (nm,
                                                             self.Param[i],
                                                             self.DObj[i],
                                                             self.DDObj[i,i],
                                                             self.Dir[i])
        return s

    def Output3(self):
        s = "FINAL VALUES\n"
        s += "   %-20s %12s\n" % ("Variable", "Value")
        s += "   %-20s %12.4e\n" % ("Obj", self.Obj)
        s += "   %-20s %12.4e\n" % ("Srel", self.Srel)
        s += "   %-20s %12.4e\n" % ("Bias", self.Bias)
        if self.HasFluct:
            s += "   %-20s %12.4e\n" % ("FluctBias", self.FluctBias)           
        for (i, nm) in enumerate(self.Names):
            if self.Fixed[i]: continue
            s += "   %-20s %12.4e\n" % (nm, self.Param[i])
        return s       

    def OutputReport(self):
        s = "========POTENTIAL SUMMARIES========\n"
        for P in self.ModSys.ForceField:
            P.SetTypeInd(None) #change the type to the default, atom type inspecific one
            s += P.ReportString()
            s += "\n\n"
        s += "========PARAMETER VALUES========\n"
        for (i, nm) in enumerate(self.Names):
            s += "   %-20s %12.4e\n" % (nm, self.Param[i])
        s += "\n"
        s += "========PARAMETER INPUT========\n"
        s += self.ModSys.ForceField.ParamString()
        s += "\n"            
        return s


           
