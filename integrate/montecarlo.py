#/usr/bin/env python


### Monte Carlo integration routine objects in SIM suite.
### coded by MSS

import numpy as np
import sim.fortran as fortran

import base.methodbase as methodbase



class MonteCarlo(methodbase.IntegrateMethodClass):
    
    CycleSource = """
>>> inputmaps
int NSteps = self.StepsPerCycle

>>> defs
int StepNum
float mcmoverannum
float Beta
float iBoxL(Dim)

>>> Main
Beta = 1.d0 / (kB * TempSet)
iBoxL = 1.d0 / BoxL
KEnergy = 0.d0
TEnergy = 0.d0

do StepNum = 0, NSteps - 1

    !draw a random number to decide which move to do
    call ran2(mcmoverannum) 
    
    [movesblock]

enddo

[external:calcenergyforces(Mode=0)]
"""

    
    ModuleVars = ["MCMoveProbs", "OldPos"]
    ModuleConsts = ["NMCMoves"]

    def __init__(self, Sys, Moves = [], StepsPerCycle = 10000):
        """Initializes a velocity verlet MD integrator."""
        object.__setattr__(self, "LibVars", {})
        methodbase.IntegrateMethodClass.__init__(self, Sys, StepsPerCycle = StepsPerCycle)
        #moves
        self.Moves = Moves
        #timestep -- doesn't mean much but in case user wants to adjust
        self.TimeStep = 1.
        #fortran subroutine for moves
        self.CycleSub = fortran.Subroutine("montecarlocycle", SubStr = "self.CycleSub", 
                                           PreCallCode = "Sys = self.Sys")
        self.CycleSub.PreLoad = self.PreLoad
        self.CycleSub.PostLoad = self.PostLoad   
        Sys.Lib.append(self.CycleSub)   

    def Reset(self):
        for Move in self.Moves:
            Move.Reset()

    def Cleanup(self):
        methodbase.IntegrateMethodClass.Cleanup(self)
        for Move in self.Moves:
            Move.Cleanup()
  
    
    def PreLoad(self):
        """Run before fortran compilation.""" 
        self.NMCMoves = len(self.Moves)
        self.MCMoveProbs = np.zeros(self.NMCMoves, dtype=float)
        self.OldPos = self.Sys.Pos.copy()
        
        #fortran code object
        fc = fortran.FortCode(Source = self.CycleSource)
        
        #initialize the moves
        s = ""
        for (i, Move) in enumerate(self.Moves):
            #let the move do preprocessing
            Move.FortPrefix = "MC%d_" % i
            Move.PreLoad()
            #add index to fortcode object
            Move.FortCode.VarIndex = i
            #add any variables to main code
            fc.AddFortCodeVars(Move.FortCode)
            #get the source code
            Source = Move.FortCode.Process()
            if i == 0:
                s += fortran.I("if (mcmoverannum <= MCMoveProbs(%d)) then" % i)
            else:
                s += fortran.I("elseif (mcmoverannum <= MCMoveProbs(%d)) then" % i)
            s += "\n"
            s += fortran.I(Source, Indent = 4)
            s += "\n"
        
        s += fortran.I("end if")
        fc.Add(s, Block = "movesblock")
        
        #update subroutine
        self.CycleSub.AddFortCode(fc)
        
        #now initialize the main variables            
        self.Sys.Lib.VarPreLoad(self, self.ModuleVars, self.ModuleConsts)
        
        
    def PostLoad(self):
        """Run after compilation."""
        self.Sys.Lib.VarPostLoad(self)        
        for Move in self.Moves:
            Move.PostLoad()
            

    def Init(self, NSteps):
        """Prepares to take MC steps."""
        self.Sys.ForceField.Eval()   
        #save the original position array
        self.OldPos = self.Sys.Pos
        #get move weights
        Weights = np.array([Move.Weight for Move in self.Moves], dtype=float)
        #get cumulative porbabilities
        self.MCMoveProbs = np.cumsum(Weights / Weights.sum())
        #move initialization
        for Move in self.Moves:
            if Move.Weight <= 0: continue
            Move.Init()
           
       
    def Cycle(self):
        """Completes one Monte Carlo cycle."""
        exec(self.CycleSub.CallObj)
        return self.StepsPerCycle

    def Finalize(self, NSteps):
        """Finalizes after MD steps."""
        for Move in self.Moves:
            if Move.Weight <= 0: continue
            Move.Finalize()
            
        
    def __setattr__(self, name, val):
        Var = self.LibVars.get(name, None)
        if Var is None:
            v = self.__dict__.get(name, None)
            if isinstance(v, np.ndarray):
                v[:] = val
            else:
                object.__setattr__(self, name, val)  
        else:
            if Var.shape:
                Var[:] = val
            else:
                Var.itemset(val)
                    
    def __getattr__(self, name):
        Var = self.LibVars.get(name, None)
        if Var is None:
            raise AttributeError(name) 
        else:
            if not Var.shape: Var = Var.item()
            return Var        


