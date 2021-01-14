#/usr/bin/env python


### Base for MD integration method objects in SIM suite.
### coded by MSS


class IntegrateMethodClass(object):

    def __init__(self, Sys, StepsPerCycle = 10):
        """Base class for integration methods."""
        self.StepsPerCycle = StepsPerCycle
        self.ErrorFlag = 0
        self.ErrorMsg = ""
        self.Sys = Sys
        
    def GetCycles(self, NSteps):
        """Converts steps to number of cycles."""
        self.CheckCycles(NSteps)
        return int(NSteps / self.StepsPerCycle)
        
    def CheckCycles(self, NSteps, OriginStr = "function"):
        """Raises an error if NSteps is not evenly divisible by StepsPerCycle."""
        if NSteps % self.StepsPerCycle > 0:
            print "WARNING: NSteps (%d) must be a multiple of StepsPerCycle (%d) in %s." % (NSteps, self.StepsPerCycle, OriginStr)        

    def Init(self, NSteps):
        """Initializes before taking integration steps."""
        #to be defined in subclasses
        pass    

    def Finalize(self, NSteps):
        """Finalizes anything after taking integration steps."""
        #to be defined in subclasses
        pass
        
    def Cycle(self):
        """Does one integration cycle; returns the number of steps."""
        #to be defined in subclasses
        pass

    def Reset(self):
        """Resets the integrator method."""
        #to be defined in subclasses
        pass

    def Cleanup(self):
        """Cleans up data structures."""
        self.Sys = None