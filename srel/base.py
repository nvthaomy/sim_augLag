#/usr/bin/env python


### Class definitions for srel routines objects in SIM suite.
### coded by MSS

import numpy as np

import sim.utility as utility
import sim.potential.base.potentialtypes as ptypes


DEBUG = False

#fractional energy difference tolerance
DiffEneFracTol = 1.e-4
#absolute energy difference tolerance in units of Sys.Units.EScale
DiffEneTol = 1.


def LogSumExp(LogTerms):
    """Computes ln(sum(exp(LogTerms))) to maintain precision."""
    LogMax = LogTerms.max()
    return LogMax + np.log(np.sum(np.exp(LogTerms - LogMax)))

def LogDiffExp(LogTerm1, LogTerm2):
    """Computes ln(sum(exp(LogTerm1)) - sum(exp(LogTerm2)) to maintain precision."""
    LogMax = max(LogTerm1, LogTerm2)
    return LogMax + np.log(np.exp(LogTerm1 - LogMax)
                           - np.exp(LogTerm2 - LogMax))

def CalcWeights(LogWeights):
    """Computes exp(LogWeights) / sum(exp(LogWeights)) to maintain precision."""
    LogMax = LogWeights.max()
    Weights = np.exp(LogWeights - LogMax)
    return Weights / Weights.sum()


def EstimateParam(Sys):
    """Estimates parameters based on arguments."""
    Sys.ForceField.EstimateGlobals()
    for P in Sys.ForceField:
        P.Estimate()
    #check bounds
    Param = Sys.ForceField.Param
    Val = Param.Val
    Min = Param.Min
    Max = Param.Max
    for i in range(len(Val)):
        if Val[i] > Max[i]:
            if Min[i] > -1.e-299:
                Val[i] = (Max[i] + Min[i]) / 2.
            else:
                Val[i] = Max[i] - 0.1 * abs(Max[i])
        elif Val[i] < Min[i]:
            if Max[i] < 1.e-299:
                Val[i] = (Max[i] + Min[i]) / 2.
            else:
                Val[i] = Min[i] + 0.1 * abs(Min[i])



def PrepSys1(Sys):
    """Prepares a system for computing energies."""
    #turn off all algorithm flags 
    Sys.Flags.CalcsOff()
    #turn off measures
    Sys.Measures.Reset()
    Sys.Measures.AllOff()
    Sys.Measures.Init()


def PrepSys2(Sys):
    """Prepares a system for computing energies and DUParam."""
    #set algorithm flags 
    Sys.Flags.CalcsOff()
    Sys.Flags.CalcDUParam = True
    #turn off measures
    Sys.Measures.Reset()
    Sys.Measures.AllOff()
    Sys.Measures.Init()


def PrepSys3(Sys):
    """Prepares a system for computing energies, DUParam, and DDUParam measures."""
    #set algorithm flags 
    Sys.Flags.CalcsOff()
    Sys.Flags.CalcDUParam = True
    #turn on DUParam measures
    Sys.Measures.Reset()
    Sys.Measures.DUParam.On()
    Sys.Measures.DDUParam.On()
    #if the step frequencies are not specified, default to every step
    if Sys.Measures.DUParam.StepFreq == 0 and Sys.Measures.DUParam.CycleFreq == 0:
        Sys.Measures.DUParam.StepFreq = 1
    if Sys.Measures.DDUParam.StepFreq == 0 and Sys.Measures.DDUParam.CycleFreq == 0:
        Sys.Measures.DDUParam.StepFreq = 1
    Sys.Measures.Init()



def CalcArgHistInt(Sys, Beta, StepsEquil, StepsProd):
    """Evaluates histograms of arguments from an integrator.
Sys : a system object with a force field defined for the model ensemble
Beta   : inverse temperature
IntStepsEqui: number of equilibration steps 
StepsProd: number of production steps 
"""
    #set temperature
    Sys.TempSet = 1./(Beta * Sys.Units.kB)
    #EQUILIBRATION PHASE
    Int = Sys.Int
    Int.Reset()
    #turn off calculation and recording of measurements
    PrepSys1(Sys)
    #run the equilibration steps
    Int.Run(StepsEquil, ProgressText = "Equilibration")
    #PRODUCTION PHASE
    Int.Reset()
    #reset argument statistics
    Sys.ForceField.Arg.ResetHist()
    #make function for evaluating  hists
    def Fn(Sys1):
        Sys1.ForceField.EvalArgHistAdd()
    #add argument distributions
    Act = Int.AddAction(Fn = Fn, StepFreq = 1)
    Int.Run(StepsProd, ProgressText = "Production")
    #delete action
    Int.DelAction(Act)
    #normalize statistics
    Sys.ForceField.Arg.NormalizeHist()
    return 
    

def GetEneTraj(Traj, Sys, ErrorNoEne = 1, ErrorDiffEne = 0, Penalties = [],
               CalcFluct = False, Beta = 0.):
    """Extracts energies stored in and evaluated from a trajectory.
Traj   : a trajectory object returning positions and potential energies
         from the target ensemble (1)
Sys : a system object with a force field defined for the model ensemble
ErrorNoEne: 1 to raise errors and 2 warnings when can't find energy in trajectory, otherwise 0
ErrorDiffEne: 1 to raise errors and 2 warnings when Traj and Sys energies differ, otherwise 0
Penalties: list of penalty objects
CalcFluct: calculate the energy fluctuations term?
Beta: 1/kBT (used if CalcFluct = True)
Returns: TrajEne, SysEne, TrajLogWeight, FluctTerm
"""   
    #save positions
    OldPos = Sys.Pos.copy()   
    
    #get the trajectory length
    n = len(Traj)

    #set all algorithm flags to none
    PrepSys1(Sys)
    if CalcFluct:
        Sys.ForceField.FluctBeta = Beta
        Sys.Flags.CalcFluct = True
        ErrorNoEne = 1
       
    #initialize the progress bar
    prog = utility.ProgressBar("Reading trajectory and evaluating energies", n)

    #go through the target trajectory and get energies    
    TrajEne = np.zeros(n, float)
    SysEne = np.zeros(n, float)
    TrajLogWeight = np.zeros(n, float)
    FluctTerm = np.zeros(n, float)
    
    #initialize any penalties
    for Pen in Penalties:
        Pen.InitializeAveraging1(n)
  
    WarnedNoEne = False
    WarnedDiffEne = False
    for (i, Pos) in enumerate(Traj):
        #get the potential energy of this configuration
        if "PEnergy" in Traj.FrameData:
            TrajEne[i] = Traj.FrameData["PEnergy"]
        elif ErrorNoEne == 2:
            if WarnedNoEne:
                print "WARNING: Could not find PEnergy in trajectory frame %d." % i
            else:
                WarnedNoEne = True
        elif ErrorNoEne == 1:
            raise KeyError("Could not find PEnergy in trajectory frame %d." % i)   
        #get the log weight of this configuration
        TrajLogWeight[i] = Traj.FrameData.get("LogWeight", 0.)
        #update the box size
        if "BoxL" in Traj.FrameData:
            Sys.ScaleBox(BoxL = Traj.FrameData["BoxL"], UpdatePos = False)
        #set the fluctuations energy
        if CalcFluct:
            Sys.ForceField.FluctE = Traj.FrameData["PEnergy"]            
        #update the model positions using the mapping function
        Sys.Pos = Pos
        #evaluate the model energy
        Sys.ForceField.Eval()
        SysEne[i] = Sys.PEnergy
        if CalcFluct:
            FluctTerm[i] = Sys.ForceField.FluctTerm
        #print out SysEne and TrajEne
        #compare energies if requested
        if ErrorDiffEne and "PEnergy" in Traj.FrameData:
            EneDiff = np.abs(TrajEne[i] - SysEne[i]) 
            FracEneDiff = EneDiff / np.abs(SysEne[i])          
            if EneDiff >= DiffEneTol * Sys.Units.EScale and FracEneDiff > DiffEneFracTol and not WarnedDiffEne:
                s = "Traj energy (%11.4e) and sys energy (%11.4e) differ in frame %d" % (TrajEne[i], SysEne[i], i)
                if ErrorDiffEne == 2:
                    print "WARNING: " + s
                    print "Suppressing further energy difference warnings."
                else:
                    raise ValueError(s)
                WarnedDiffEne = True
        EneDiff = np.abs(TrajEne[i] - SysEne[i])
        FracEneDiff = EneDiff / np.abs(SysEne[i])
 
        if ErrorDiffEne == 2 and i%100==0:
            print ('%d %11.4e %11.4e %11.4e %11.4e' %(i, TrajEne[i], SysEne[i], EneDiff, FracEneDiff))
        #update any penalties
        if len(Penalties):
            Sys.Measures.Eval()
            for Pen in Penalties:
                Pen.UpdateAveraging1(i)
        #update the progress bar
        prog.Update(i)
    prog.Clear()   

    #restore positions
    Sys.Pos = OldPos

    #turn off stuff
    PrepSys1(Sys) 
    
    return TrajEne, SysEne, TrajLogWeight, FluctTerm  


def ParseArgRanges(Traj, Sys):
    """Gets argument ranges from a trajectory.
Traj: a trajectory object 
Sys : a system object with a force field defined for the model ensemble
Returns a list of arg min and max arrays
"""
   
    #get the trajectory length
    n = len(Traj)

    #initialize the progress bar
    prog = utility.ProgressBar("Examining trajectory argument mins and maxes", n)

    #go through the target trajectory and get argument mins and maxes
    Sys.ForceField.Arg.ResetStats()
    for (i, Pos) in enumerate(Traj):
        #update the box size
        BoxL = Traj.FrameData.get("BoxL", None)
        if not BoxL is None:
            Sys.ScaleBox(BoxL = BoxL, UpdatePos = False)
        #update the model positions 
        Sys.Pos = Pos
        #evaluate the model args
        Sys.ForceField.EvalArgStatsAdd()
        #update the progress bar
        prog.Update(i)  
    prog.Clear()

    #normalize statistics
    Sys.ForceField.Arg.NormalizeStats()
    
    if DEBUG:
        print "Argument parameters min / max / histnbins:"
        for P in Sys.ForceField:
            print "Force field component %s" % P.Name
            print " ", P.Arg.Min, P.Arg.Max, len(P.Arg.Hist) 
            
    return
    

def ParseArgHist(Traj, Sys, Beta, Reweight = False):
    """Parses argument distributions from a trajectory.
Traj: a trajectory object 
Sys : a system object with a force field defined for the model ensemble
Beta   : inverse temperature
Reweight: reweight trajectory using differences between Sys and Traj enes
Returns a representative position array and box size from the system from the trajectory.
"""
   
    #get the trajectory length
    n = len(Traj)

    #save positions
    OldPos = Sys.Pos.copy()  
    
    #get the trajectory and system energies
    ret = GetEneTraj(Traj, Sys, ErrorDiffEne = 0, ErrorNoEne = Reweight)   
    TrajEne, SysEne, TrajLogWeight, FluctTerm = ret
    
    #prepping for reweighting    
    if Reweight:
           
        #compute model weights
        LogWeight = TrajLogWeight + Beta * (TrajEne - SysEne)
        Weights = CalcWeights(LogWeight)
        
    else:
          
        #compute weights for trajectory
        Weights = CalcWeights(TrajLogWeight)  
         
    #find a reference configuration
    RefInd = None
    if RefInd is None:
        WeightAvg = np.mean(Weights)
        WeightStd = np.std(Weights)
        if np.all(WeightStd < 1.e-6 * WeightAvg):
            #weights are all about the same
            #just choose the last configuration
            RefInd = n - 1
        else:
            #choose a configuration with the maximum weight
            RefInd = np.argmax(Weights)

    #initialize the progress bar
    prog = utility.ProgressBar("Making argument histograms", n)

    #go through the target trajectory and get argument distributions
    Sys.ForceField.Arg.ResetHist()
    RefPos = None
    RefBoxL = None
    for (i, Pos) in enumerate(Traj):        
        #get the box size
        BoxL = Traj.FrameData.get("BoxL", None)
        #update the box size
        if not BoxL is None:
            Sys.ScaleBox(BoxL = BoxL, UpdatePos = False)
        #update the model positions 
        Sys.Pos = Pos
        #evaluate the model hist
        Sys.ForceField.EvalArgHistAdd(Weight = Weights[i])
        #find the reference configuration
        if i == RefInd:
            RefPos = Pos.copy()
            if not BoxL is None:
                RefBoxL = BoxL.copy()
        #update the progress bar
        prog.Update(i)
        
    prog.Clear()

    #normalize the histogram per frame
    Sys.ForceField.Arg.NormalizeHist()

    Sys.Pos = OldPos   

    return RefPos, RefBoxL    
    


def CalcAvgsTraj(Traj, Sys, Beta, Reweight = False, CalcDeriv = True, 
                 Penalties = [], CalcFluct = False):
    """Computes averages from a trajectory.
Traj   : a trajectory object returning positions and potential energies
         from the target ensemble (1)
Sys : a system object with a force field defined for the model ensemble
Beta   : inverse temperature
Reweight: True to reweight based on energy differences between the
          trajectory and computed with Sys
CalcDeriv: True to compute Srel derivatives; otherwise zero arrays returned
Penalties: list of penalty objects
CalcFluct: add in the fluctuation calculations
Returns: DUParam[:], SqDUParam[:,:], DDUParam[:,:], Count[:], ReweightFrac, SrelReweighTerm, AvgPEnergy, AvgFluctTerm
"""

    #save positions
    OldPos = Sys.Pos.copy()    

    #get the trajectory length
    n = len(Traj)
    
    #get parameters
    N = Sys.ForceField.Param.N
    
    #tell penalties how to update
    for Pen in Penalties:
        Pen.CalcDeriv = CalcDeriv
    
    #get the trajectory and system energies
    ret = GetEneTraj(Traj, Sys, ErrorDiffEne = 0, ErrorNoEne = 1, Penalties = Penalties,
                     CalcFluct = CalcFluct, Beta = Beta)
    TrajEne, SysEne, TrajLogWeight, FluctTerm = ret
  
    if Reweight:    
        
        #compute model weights without trajectory to get reweight fraction
        LogWeight = Beta * (TrajEne - SysEne)
        Weight = CalcWeights(LogWeight)
        ReweightFrac = np.exp(-np.sum(np.where(Weight > 0., Weight * np.log(Weight), 0.))) / len(Traj)
        
        #compute full model weights
        LogWeight = TrajLogWeight + Beta * (TrajEne - SysEne)
        Weight = CalcWeights(LogWeight)
        
        #compute srel terms
        Delta = Beta * (TrajEne - SysEne)
        LogNum = LogSumExp(Delta + TrajLogWeight)
        LogDem = LogSumExp(TrajLogWeight)
        SrelReweighTerm = LogNum - LogDem

    else:
               
        #compute weights for trajectory
        Weight = CalcWeights(TrajLogWeight)  
        
        #compute reweighting fraction
        ReweightFrac = 1.
        
        #compute Srel terms
        SrelReweighTerm = 0.
        
    #compute average energy
    AvgPEnergy = np.sum(SysEne * Weight)
    #compute average fluctuation term
    AvgFluctTerm = np.sum(FluctTerm * Weight)
    
    #finalize averaging for any penalties
    for Pen in Penalties:
        Pen.FinalizeAveraging1(Weight)
    
    #compute effective number of counts
    #uses error propagation rules for "standard error of a weighted mean"
    Count = np.ones((N,), dtype = float) / np.sum(Weight * Weight)
    
    #check if we don't want derivatives
    if not CalcDeriv:
        #restore positions
        Sys.Pos = OldPos
        #turn off stuff
        PrepSys1(Sys) 
        #create return variables
        DUParam = np.zeros((N,), float)
        SqDUParam = np.zeros((N,N), float)
        DDUParam = np.zeros((N,N), float)
        return DUParam, SqDUParam, DDUParam, Count, ReweightFrac, SrelReweighTerm, AvgPEnergy, AvgFluctTerm
        
    #initialize progress bar for second pass
    prog = utility.ProgressBar("Calculating trajectory averages", n)

    #turn on parameter derivative calculation
    PrepSys3(Sys)
    if CalcFluct:
        Sys.ForceField.FluctBeta = Beta
        Sys.Flags.CalcFluct = True 
    
    #initialize penalties
    for Pen in Penalties:
        Pen.InitializeAveraging2()       
    
    #go through the target trajectory again and average derivatives
    for (i, TrajPos) in enumerate(Traj):
        #update the box size
        if "BoxL" in Traj.FrameData:
            Sys.ScaleBox(BoxL = Traj.FrameData["BoxL"], UpdatePos = False)
        #set the fluctuations energy
        if CalcFluct:
            Sys.ForceField.FluctE = Traj.FrameData["PEnergy"]             
        #update the model positions
        Sys.Pos = TrajPos
        #evaluate the model energy
        Sys.ForceField.Eval()
        #update the measures
        Sys.Measures.Eval(Weight = Weight[i])
        #update penalties
        for Pen in Penalties:
            Pen.UpdateAveraging2(Weight[i])
        #update the progress bar
        prog.Update(i)
    prog.Clear()

    #finalize penalties
    for Pen in Penalties:
        Pen.FinalizeAveraging2(Beta)

    #get parameter averages
    DUParam = Sys.Measures.DUParam.Avg
    SqDUParam = Sys.Measures.DUParam.AvgSq 
    DDUParam = Sys.Measures.DDUParam.MatrixAvg

    #set the derivative for fixed parameters to zero    
    Fixed = Sys.ForceField.Param.Fixed
    DUParam[Fixed] = 0.
    SqDUParam[np.ix_(Fixed,Fixed)] = 0.  
    DDUParam[np.ix_(Fixed,Fixed)] = 0.    

    #restore positions
    Sys.Pos = OldPos  

    #turn off stuff
    PrepSys1(Sys)  
   
    return DUParam, SqDUParam, DDUParam, Count, ReweightFrac, SrelReweighTerm, AvgPEnergy, AvgFluctTerm    


def CalcAvgsHist(Sys, CalcDeriv = True):
    """Evaluates parameter derivatives for a trajectory.  Each call
uses pre-calculated histograms.
Sys : a system object with a force field defined for the model ensemble
CalcDeriv: True to compute Srel derivatives; otherwise zero arrays returned
Returns DUParam, DDUParam, AvgPEnergy
"""
    #check that we can do it with this potential
    for P in Sys.ForceField:
        if P.NoHistEval:
            raise TypeError("Cannot evaluate energies with histograms for %s." % P.Name)
        elif P.Arg.HistNBin == 0:
            raise TypeError("Histogram for %s has zero bins." % P.Name)
    
    #turn on parameter calculation
    if CalcDeriv:
        PrepSys2(Sys)
    else:
        PrepSys1(Sys)

    #evaluate using the histograms
    Sys.ForceField.EvalArg()
    AvgPEnergy = Sys.PEnergy

    #get the derivatives
    #also set the derivative for fixed parameters to zero   
    if CalcDeriv:
        DUParam = Sys.ForceField.Param.DU.copy()
        DDUParam = Sys.ForceField.Param.DDU.Matrix 
        Fixed = Sys.ForceField.Param.Fixed
        DUParam[Fixed] = 0.
        DDUParam[np.ix_(Fixed,Fixed)] = 0.
    else:
        N = Sys.ForceField.Param.N
        DUParam = np.zeros(N, dtype=float)
        DDUParam = np.zeros((N,N), dtype=float)

    Sys.ForceField.Param.Fixed
    
    #turn stuff off
    PrepSys1(Sys)

    return DUParam, DDUParam, AvgPEnergy


def CalcSrelTraj(TarTraj, Sys, Beta, CalcDeriv = True, Penalties = [], CalcFluct = False):
    """Evaluates the relative entropy and derivatives entirely from a target
trajectory.  Each call requires two sorts through the trajectory.
TarTraj: a mapped trajectory object returning positions from the model
         ensemble and and potential energies from the target
Sys : a system object with a force field defined for the model ensemble
Beta   : inverse temperature
CalcDeriv : True to compute derivatives
Penalties: list of penalty objects
CalcFluct: add in the fluctuation calculations
Returns: Srel, DSrel[:], DDSrel[:,:], AvgFluctTerm
"""

    #save positions
    OldPos = Sys.Pos.copy()   
    
    #get the trajectory length
    n = len(TarTraj)
    
    #get parameters
    N = Sys.ForceField.Param.N
    N2 = Sys.ForceField.Param.N2
    
    #tell penalties how to update
    for Pen in Penalties:
        Pen.CalcDeriv = CalcDeriv
    
    #get the energies in the trajectory    
    ret = GetEneTraj(TarTraj, Sys, ErrorNoEne = 1, ErrorDiffEne = 0, Penalties = Penalties,
                     CalcFluct = CalcFluct, Beta = Beta)  
    TarEne, ModEne, TarLogWeight, FluctTerm = ret
   
    #compute target weights
    TarWeight = CalcWeights(TarLogWeight)

    #compute the relative entropy
    Delta = Beta * (TarEne - ModEne)
    LogNum = LogSumExp(Delta + TarLogWeight)
    LogDem = LogSumExp(TarLogWeight)
    Srel = LogNum - LogDem
    
    #compute model weights   
    ModLogWeight = Delta + TarLogWeight - Srel
    ModWeight = CalcWeights(ModLogWeight)
    
    #compute the average fluctuation
    AvgFluctTerm = np.sum(TarWeight * FluctTerm)
    
    #finalize averaging for any penalties
    for Pen in Penalties:
        Pen.FinalizeAveraging1(ModWeight)

    #check if all we wanted was to compute Srel
    if not CalcDeriv:
        Sys.Pos = OldPos
        return Srel, np.zeros((N,), float), np.zeros((N,N), float)

    #turn on parameter derivative calculation
    PrepSys2(Sys)
    if CalcFluct:
        Sys.ForceField.FluctBeta = Beta
        Sys.Flags.CalcFluct = True
      
    #initialize progress bar for second pass
    prog = utility.ProgressBar("Calculating trajectory averages", n)

    #go through the target trajectory again and average derivatives
    ModDUParam = np.zeros(N, float)
    ModSqDUParam = np.zeros((N,N), float)
    ModDDUParam = np.zeros(N2, float)
    TarDUParam = np.zeros(N, float)
    TarDDUParam = np.zeros(N2, float)
    Param = Sys.ForceField.Param
    
    #initialize penalties
    for Pen in Penalties:
        Pen.InitializeAveraging2()

    for (i, Pos) in enumerate(TarTraj):
        #update the box size
        if "BoxL" in TarTraj.FrameData:
            Sys.ScaleBox(BoxL = TarTraj.FrameData["BoxL"], UpdatePos = False)         
        #update the model positions 
        Sys.Pos = Pos
        #evaluate the model energy
        Sys.ForceField.Eval()
        #update the running averages for each parameter
        ModDUParam += ModWeight[i] * Param.DU
        ModSqDUParam += ModWeight[i] * np.outer(Param.DU, Param.DU)
        ModDDUParam += ModWeight[i] * Param.DDU.Array
        TarDUParam += TarWeight[i] * Param.DU
        TarDDUParam += TarWeight[i] * Param.DDU.Array
        #update penalties
        for Pen in Penalties:
            Pen.UpdateAveraging2(ModWeight[i])
        #update the progress bar
        prog.Update(i)
    prog.Clear()
    
    #finalize penalties
    for Pen in Penalties:
        Pen.FinalizeAveraging2(Beta)

    #convert to 2D matrices 
    ModDDUParam = Sys.ForceField.Param.DDU.MakeMatrix(ModDDUParam)
    TarDDUParam = Sys.ForceField.Param.DDU.MakeMatrix(TarDDUParam)

    #compute the Srel derivatives
    DSrel = Beta * (TarDUParam - ModDUParam)
    DDSrel = Beta * (TarDDUParam - ModDDUParam + Beta * ModSqDUParam 
                     - Beta * np.outer(ModDUParam, ModDUParam))
    #set the derivative for fixed parameters to zero
    Fixed = Sys.ForceField.Param.Fixed
    if not Fixed is None:
        DSrel[Fixed] = 0.
        DDSrel[np.ix_(Fixed,Fixed)] = 0.

    Sys.Pos = OldPos

    #turn stuff off
    PrepSys1(Sys)

    return Srel, DSrel, DDSrel, AvgFluctTerm    

        
    
def CheckAllModSys(AllModSys):
    """Checks to make sure that all systems have same parameters."""
    ff0 = AllModSys[0].ForceField
    for ModSys in AllModSys[1:]:
        ff = ModSys.ForceField
        if not len(ff0) == len(ff):
            raise ValueError("Systems' force fields are not the same.")
        for (ff0i, ffi) in zip(ff0, ff):
            if not type(ff0i) == type(ffi):
                raise ValueError("Systems' force fields are not the same.")
        if not len(ff0.Param.Val) == len(ff.Param.Val):
            raise ValueError("Systems' force fields are not the same.")
        if not np.all(ff0.Param.Fixed == ff.Param.Fixed):
            raise ValueError("Fixed parameters must be the same between systems.")
            
    
    
    
    
        
            
        
        
    
            
    
