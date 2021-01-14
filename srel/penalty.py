#/usr/bin/env python


### Penalty constraints for relative entropy optimizer.
### coded by MSS


import numpy as np


class PenaltyClass(object):
    
    def __init__(self, Sys, Measure, Target, Coef = 1., Name = "penalty", 
                 MeasureScale = 1., ValInd = 0):
        """Initializes a penalty constraint in the minimization, of the form
Coef * (<Measure> - Target)^2  where <Target> is the ensemble average from Measure.  
Coef can be gradually increased over successive minimizations until the constraint
is satisfied within reasonable tolerance.
Sys: The system of interest
Measure:  An instance of a MeasureAvgClass object that is used to compute 
          the average to be constrained 
Target:  The target value of the average
Coef:  The coefficient of the penalty function when added to the main objective
Name:  An optional name
ValInd: Index of the value to use in measure, if more than one
MeasureScale: optional scale factor for the measurement and target value"""
        if Measure == "PEnergy" or getattr(Measure, "Name", "") == "PEnergy":
            self.__Mode = 1
            self.Measure = Sys.Measures.PEnergy
        elif Measure == "Virial" or getattr(Measure, "Name", "") == "Virial":
            self.__Mode = 2
            self.Measure = Sys.Measures.Virial
        else:
            self.__Mode = 0
            from sim.measure.base import MeasureAvgClass
            if not isinstance(Measure, MeasureAvgClass):
                raise TypeError("MeasureAvgObj must be an instance of MeasureAvgClass or subclass.")
            if not Measure.Sys is Sys:
                raise ValueError("Measure does not appear to be associated with Sys.")
            if ValInd >= Measure.NVal:
                raise ValueError("ValInd is greater than number of values in Measure.")
            self.Measure = Measure
        self.Target = Target
        self.Coef = Coef
        self.LagMult = 0.
        self.Name = Name
        self.ValInd = ValInd
        self.Obj = 0.
        self.MeasureScale = MeasureScale
        self.KeepDerivs = False
        self.CalcDeriv = False
        
    def InitializeOptimization(self):
        self.LagMult = 0.
        self.Obj = 0.
        
    def InitializeAveraging1(self, n):
        if self.CalcDeriv: return
        self.Measure.On()
        self.Measure.ResetAvgs()
        self.ValArray = np.zeros(n, dtype=float)
        
    def UpdateAveraging1(self, i):
        if self.CalcDeriv: return
        if self.Measure.NVal == 1:
            self.ValArray[i] = self.Measure.Val
        else:
            self.ValArray[i] = self.Measure.Val[self.ValInd]
            
    def FinalizeAveraging1(self, Weight):
        if self.CalcDeriv: return
        self.DAvg = None
        self.DDAvg = None
        self.Avg = np.sum(self.ValArray * Weight)
        del self.ValArray
    
    def InitializeAveraging2(self):
        if not self.CalcDeriv: return
        self.Measure.On()
        self.Measure.ResetAvgs()
        NParam = len(self.Measure.Sys.ForceField.Param.Val)
        self.WeightSum = 0.
        self.DUParamA = np.zeros((NParam), dtype=float)
        self.DUParamSqV = np.zeros((NParam, NParam), dtype=float)
        self.DDUParamA = np.zeros((NParam, NParam), dtype=float)
        if self.__Mode == 2:
            Sys = self.Measure.Sys
            Sys.Measures.DWParam.On()
            Sys.Measures.DDWParam.On()
            Sys.Measures.DUParamDWParam.On()
        
    def UpdateAveraging2(self, Weight):
        if not self.CalcDeriv: return
        if self.Measure.NVal == 1:
            Val = self.Measure.Val
        else:
            Val = self.Measure.Val[self.ValInd]
        Sys = self.Measure.Sys
        DUParam = Sys.ForceField.Param.DU
        DDUParam = Sys.ForceField.Param.DDU.Matrix
        self.DUParamA += (Val * DUParam * Weight)
        self.DUParamSqV += (Val * np.outer(DUParam, DUParam) * Weight)
        self.DDUParamA += (Val * DDUParam * Weight)
        self.WeightSum += Weight
        
    def FinalizeAveraging2(self, Beta):
        if not self.CalcDeriv: return   
        Sys = self.Measure.Sys
        self.Avg = self.Measure.Avg
        Avg = self.Avg
        w = self.WeightSum + 1.e-300
        DUParamA = Beta * self.DUParamA / w
        del self.DUParamA
        DUParamSqV = Beta * Beta * self.DUParamSqV / w
        del self.DUParamSqV
        DDUParamA = Beta * self.DDUParamA / w
        del self.DDUParamA
        DUParam = Beta * Sys.Measures.DUParam.Avg
        DUParamSq = Beta * Beta * Sys.Measures.DUParam.AvgSq
        DDUParam = Beta * self.Measure.Sys.Measures.DDUParam.MatrixAvg
        self.DAvg = DUParam * Avg - DUParamA
        self.DDAvg = (2. * np.outer(DUParam, DUParam) * Avg
                      - np.outer(DUParam, DUParamA) - np.outer(DUParamA, DUParam)
                      - DUParamSq * Avg + DUParamSqV
                      + DDUParam * Avg - DDUParamA)
        if self.__Mode == 1:
            self.DAvg += DUParam / Beta
            self.DDAvg += (2 * (np.outer(DUParam, DUParam) - DUParamSq) + DDUParam) / Beta
        elif self.__Mode == 2:
            DWParam = Sys.Measures.DWParam.Avg
            DDWParam = Sys.Measures.DDWParam.MatrixAvg
            DUParamDWParam = Beta * Sys.Measures.DUParamDWParam.MatrixAvg
            self.DAvg += DWParam
            a = np.outer(DWParam, DUParam) - DUParamDWParam
            self.DDAvg += a + a.T + DDWParam
        
    def UpdateObj(self, Obj, Bias, DObj, DDObj):
        """Updates the contributions to Obj, Bias, DObj, and DDObj."""
        Delta = self.MeasureScale * (self.Avg - self.Target)
        print('New Calculate Avg value: {}'.format(self.Avg))
        print('Target: {}'.format(self.Target))
        print('Delta: {}'.format(Delta))
        self.Obj = 0.5 * self.Coef * np.sum(Delta**2) - self.LagMult * Delta
        print('Obj {}'.format(self.Obj))
        Obj += self.Obj
        Bias += self.Obj
        if not self.DAvg is None:
            ScDAvg = self.MeasureScale * self.DAvg
            ScDDAvg = self.MeasureScale * self.DDAvg
            DObj += self.Coef * Delta * ScDAvg  - self.LagMult * ScDAvg
            DDObj += self.Coef * (Delta * ScDDAvg + np.outer(ScDAvg, ScDAvg)) - self.LagMult * ScDDAvg
            if not self.KeepDerivs:
                self.DAvg = None
                self.DDAvg = None
        return Obj, Bias, DObj, DDObj
        
    def UpdateLagMult(self):
        self.LagMult = self.LagMult - self.Coef * self.MeasureScale * (self.Avg - self.Target)
        
    def CleanUp(self):
        self.Measure = None
        

