#/usr/bin/env python


### Base for MD integration routine objects in SIM suite.
### coded by MSS

import methodbase

import sim.utility as utility

import time
import numpy as np


class Action(object):
    """Base class for actions."""

    def __init__(self, StepFreq = 0, CycleFreq = 0, TimeFreq = 0, 
                 Fn = None, InitFn = None, FinalFn = None, Name = None):
        """Initializes an action for integrations."""
        self.StepFreq = StepFreq
        self.CycleFreq = CycleFreq
        self.TimeFreq = TimeFreq
        self.Name = Name
        def Do(Sys):
            """Runs the action."""
            #to be defined in subclasses
            pass
        if Fn is None:
            self.Do = Do
        else:
            self.Do = Fn
        self.InitFn = InitFn
        self.FinalFn = FinalFn
        self.Reset()
        
    def Reset(self):
        """Resets counters in action."""
        self.LastStepNum = 0
        self.LastTime = time.time()
        
    def __repr__(self):
        """Returns Action name"""    
        if self.Name:
            return "<%s>" % self.Name
        else:
            return repr(self.__class__)


class MethodsContainer(list):
    def __init__(self):
        pass

    
class Integrator(object):

    def __init__(self, Sys, ProgressText = None):
        """Base class for integrators."""
        self.Sys = Sys
        self.Actions = []
        self.Methods = MethodsContainer()
        self.Method = None
        self.ProgressText = ProgressText
        self.Reset()

    def Reset(self):
        "Resets counters and method."
        self.StepNum = 0
        self.CycleNum = 0
        if not self.Method is None:
            self.Method.Reset()
        for Act in self.Actions:
            Act.Reset()

    def AddMethod(self, MethodClass, *args, **kwargs):
        """Adds an integration method."""
        if not issubclass(MethodClass, methodbase.IntegrateMethodClass):
            raise TypeError("MethodClass not an IntegrateMethodClass type.")
        MethodObj = MethodClass(self.Sys, *args, **kwargs)
        self.Methods.append(MethodObj)
        setattr(self.Methods, MethodClass.__name__, MethodObj)
        self.Method = MethodObj
        return MethodObj

    def AddAction(self, Act = None, **kwargs):
        if isinstance(Act, Action):
            self.Actions.append(Act)
        else:
            Act = Action(**kwargs)
            self.Actions.append(Act)
        return Act

    def DelAction(self, Act):
        self.Actions = [a for a in self.Actions if not a is Act]
        
    def Cleanup(self):
        self.Sys = None
        for (name,var) in self.Methods.__dict__.items():
            if isinstance(var, methodbase.IntegrateMethodClass):
                var.Cleanup()           

    def DoAllActions(self):
        """Runs all actions."""
        for Act in self.Actions:
            Act.Do(self.Sys)

    def Run(self, NSteps, ProgressText = None, Reset = False):
        """Runs for a preset number of integration steps."""
        #record start time
        self.StartTime = time.time()
        #check the system
        self.Sys.Check()
        #initialize
        self.Method.Init(NSteps)
        #get number of cycles
        NCycles = self.Method.GetCycles(NSteps)
        #setup progress bar
        if ProgressText is None:
            ProgressText = self.ProgressText
        if not ProgressText is None:
            PB = utility.ProgressBar(ProgressText, Steps = NCycles)        
        #initialize measures
        self.Sys.Measures.Init()
        for Measure in self.Sys.Measures:
            if Measure.Active:
                self.Method.CheckCycles(Measure.StepFreq, Measure.Name)
        #initialize actions
        for Act in self.Actions:
                if not Act.InitFn is None and (Act.StepFreq > 0 or Act.CycleFreq > 0 or Act.TimeFreq > 0):
                    Act.InitFn(self.Sys)
                    self.Method.CheckCycles(Act.StepFreq, str(Act))                    
        #run the cycles
        for i in xrange(NCycles):
            #run one cycle
            StepsTaken = self.Method.Cycle()
            #check errors
            if self.Method.ErrorFlag:
                print "Error found during integration step.  Msg is:\n" + self.Method.ErrorMsg
                raise StandardError("Integration error")
            #update number of steps
            self.StepNum += StepsTaken
            self.CycleNum += 1
            #evaluate the measures
            self.Sys.Measures.Eval(self.StepNum, self.CycleNum)
            #do actions
            Time = time.time()
            for Act in self.Actions:
                if (Act.StepFreq > 0 and self.StepNum - Act.LastStepNum >= Act.StepFreq) \
                   or (Act.CycleFreq > 0 and self.CycleNum % Act.CycleFreq == 0) \
                   or (Act.TimeFreq > 0 and Time - Act.LastTime > Act.TimeFreq):
                    Act.Do(self.Sys)
                    Act.LastTime = Time
                    Act.LastStepNum = self.StepNum
            if not ProgressText is None:
                PB.Update(i)
        #finalize method
        self.Method.Finalize(NSteps)
        #finalize measures
        self.Sys.Measures.Finalize()
        #finalize actions
        for Act in self.Actions:
                if not Act.FinalFn is None and (Act.StepFreq > 0 or Act.CycleFreq > 0 or Act.TimeFreq > 0):
                    Act.FinalFn(self.Sys)
        if not ProgressText is None:
            PB.Clear()
        #check for NaN errors
        if np.any(np.isnan(self.Sys.Pos)):
            print "Error found after integration.  NaN found in Sys.Pos."
            raise StandardError("Integration error.")
        if np.any(np.isnan(self.Sys.Vel)):
            print "Error found after integration.  NaN found in Sys.Vel."
            raise StandardError("Integration error.")
        #if we are to reset cycle numbers, step numbers, etc. do that now
        if Reset:
            self.Reset()
        #record stop time
        self.StopTime = time.time()
        #calculate duration
        self.TimeElapsed = self.StopTime - self.StartTime
    
        



            
            


            
    
        
            
        
        
    
            
    