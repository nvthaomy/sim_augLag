#/usr/bin/env python


### Relative entropy optimizer for saved trajectories.
### coded by MSS


import numpy as np
import time
import atexit
import os
import pickle

import sim.traj as traj
import sim.utility as utility

import optimize
import base


#debug args
DEBUGHESSIAN = False

#format for plots
EnablePlots = True
PlotFmt = "pdf"


#search direction modes
DirModes = {"SD" : "STEEPEST DESCENT",
            "CG"  : "CONJUGATE GRADIENT",
            "HD"  : "HESSIAN DESCENT"}

#max values for dx for each mode
dxMaxMode = {"SD" : 1.e300,
             "CG"  : 1.e300,
             "HD"  : 1.0}

#initial guesses for dx
dxModeDefaults = {"SD" : 0.01,
                  "CG"  : 0.10,
                  "HD"  : 0.10}
                        
#Staged tolerances for Augmented Larangian method
StageCGFracTols = None
StageCGGradTols = None

#Default options and parameters for minimization procedures
Defaults = {"dxScaleUp" : 1.5,
            "dxScaleDown" : 0.2,
            "MinReweightFrac" : 0.15,
            "MinReweightFrames" : None,
            "MinDeltaReweightFrac" : 0.05,
            "MaxDeltaReweightFrac" : 0.10,
            "ReweightTar" : False,  # False to reweight target trajectory directly without ever making
                                    # a model trajectory.  Rarely use, only if very slight CGing. 
            "Autodx": True,
            "AutodxErrMin" : 1.0,
            "AutodxFracErrMin" : 0.05,
            "AutodxErrMax" : 2.0,
            "AutodxFracErrMax" : 0.1,
            "LSFracTol" : 1.e-7,
            "LSAbsTol" : 1.e-7,
            "LSAccel" : 2.0,
            "LSSlow" : 0.75, 
            "LSBackScale" : 0.1,
            "LSMinBackScale" : 1.e-8,
            "LSMaxIncr" : 100.,
            "LSNumTol" : 1.e-8,   
            "LSVerbose" : False, 
            "LSMaxIter" : 200,
            "CGFracTol" : 1.e-6,
            "CGAbsTol" : 1.e-6,
            "CGGradTol" : 1.e-2,
            "CGFormula" : 1,
            "CGRestartIters" : 20,
            "CGRestartCount" : 0,
            "CGBeta" : 0., 
            "NRAbsTol" : 1.e-8,
            "NRFracTol" : 1.e-8,
            "NRGradTol" : 1.e-1,
            "FilterHessian" : False,
            "HessianMaxDir" : 1.e20,
            "HessianPad" : 1.e-5,
            "UseHessian" : True,
            "MinEig" : 0.,
            "ReweightOldModTraj" : True,
            "MaxFracChange" : 0.5,
           }


def MakeModTraj(ModSys, FilePrefix, StepsEquil, StepsProd, ElecSys, StepsStride = 1, Verbose = False):
    """Makes a trajectory for the model system and returns ModTraj.
Note: ModTraj must expose a FileName attribute and the Delete() method."""   
    #equilibrate
    ModSys.Int.Run(StepsEquil, ProgressText = "Equilibration")
    #open trajectory for writing
    ModTrajFn = FilePrefix + ".traj.dat"
    ModTraj = traj.RawWrite(ModTrajFn)
    #add action to save to trajectory
    ModTraj.AddAction(ModSys.Int, StepFreq = StepsStride)
    #run production
    ModSys.Int.Run(StepsProd, ProgressText = "Production")
    #remove trajectory saving action
    ModTraj.DelAction(ModSys.Int)
    #close trajectory
    ModTraj.Close()   
    return traj.Raw(ModTrajFn) 
    
    

class OptimizeTrajClass(optimize.OptimizeClass):
    
    def __init__(self, ModSys, Map, Traj, Beta = None,
                 FilePrefix = None,
                 Verbose = True,
                 TempFileDir = None,
                 UseTarHists = True,
                 UseTrajBoxL = True,
                 RefPosInd = None,
                 SaveLoadArgData = True, 
                 BlurSigma = None,
                 BlurNSample = 1, 
                 ElecSys = None,        #MN 2020.01.05       
                 **kwargs):
        """Initializes a Srel-minimization class based on a reference target trajectory.
ModSys: System object for the model (coarse-grained) system
Map: PosMap object that maps target (atomic) to model (CG) coordinates
Traj: trajectory object returning positions in target system; may be a list
Beta: inverse temperature * kB; defaults to that of ModSys
BoxL: box lengths of system; defaults to that of ModSys
FilePrefix: file prefix for recording iteration history, histograms, etc
Verbose: flag for detailed output; defaults to True
UseTarHists: flag for using histograms of potential energy arguments for
             target avgs rather than explicit traj reprocessing
UseTrajBoxL: True to use reference box size from trajectory
RefPosInd: integer index of reference configuration for initial MD runs; None to auto-select
SaveLoadArgData: True to save/load argument data for fast reloading
BlurSigma: degree of blurring of atomistic trajectory
BlurNSample: number of blur frames to sample for each atomistic frame
"""
        #MN 2020.01.05
        self.ElecSys = ElecSys

        #set temperature and box size
        if Beta is None:
            self.Beta = 1. / (ModSys.Units.kB * ModSys.TempSet)
        else:
            self.Beta = Beta
            ModSys.TempSet = 1. / (ModSys.Units.kB * Beta)
        Pos = Traj[0]
        TrajBoxL = Traj.FrameData.get("BoxL", None)
        if UseTrajBoxL:
            if TrajBoxL is None:
                print "WARNING: No box size in trajectory to set system size; using current in system object."
            else:
                ModSys.ScaleBox(TrajBoxL)   
        elif not TrajBoxL is None:
            if not np.allclose(TrajBoxL, ModSys.BoxL):
                print "WARNING: Trajectory has box size that differs from current system box size."
            
        #initialize class
        optimize.OptimizeClass.__init__(self, ModSys, ElecSys,**kwargs)
        
        #save system settings
        self.TempSet = ModSys.TempSet
        self.RefBoxL = ModSys.BoxL.copy()
        #set file names
        self.FilePrefix = FilePrefix
        self.TempFilePrefix = None
        self.TempFileDir = TempFileDir
        
        #set verbose options
        self.Verbose = Verbose
        #set tar hists option
        if self.HasFluct and UseTarHists:
            raise ValueError("Cannot have UseTarHists=True when system has a fluctuation potential.")
        self.UseTarHists = UseTarHists
        #save and load target argument data 
        self.SaveLoadArgData = SaveLoadArgData
        #add default options to class
        for (attr, val) in Defaults.items():
            setattr(self, attr, val)
                   
        #reference position index
        self.RefPosInd = RefPosInd
        self.RefPos = None
        self.UseTrajBoxL = UseTrajBoxL
        
        #saved argument data (to be computed)
        self.SavedTarArgData = None
        self.TarIsParsed = False
        
        #map the trajectory atomic coords into the CG coords
        MappedTraj = traj.Mapped(Traj, Map, BoxL = self.RefBoxL)
        #convert the trajectory, keeping weight and energy frame data
        if self.UseTarHists and BlurSigma is None:
            #don't convert trajectory because we will just always use histograms anyways
            self.TarTrajFn = ""
            self.TarTraj = MappedTraj
        else:
            #preprocess the trajectory into "raw" format for faster access
            self.TarTrajFn = utility.NewTempFile(Suffix=".dat", Dir=self.TempFileDir)
            if self.Verbose:
                print "Preprocessing target trajectory into file %s"  % self.TarTrajFn
            traj.Convert(MappedTraj, traj.RawWrite, self.TarTrajFn, Verbose = Verbose,
                         FrameDataFields = ["LogWeight", "PEnergy"],
                         BlurSigma = BlurSigma, BlurNSample = BlurNSample)
            self.TarTraj = traj.Raw(self.TarTrajFn)
            
        #set the model traj
        self.ModTraj = None
        self.OldModTraj = None
        self.NewModTraj = None
        
        #set the trajectory-making function
        self.MakeModTrajFunc = MakeModTraj
        
        #automatically delete the preprocessed trajectory at exit
        def AtExitFunc():
            self.Cleanup()
        atexit.register(AtExitFunc)
                
        #check that histograms have been defined   
        if self.UseTarHists:
            for P in self.ModSys.ForceField:
                if P.Arg.HistNBin == 0:
                    raise ValueError("Must define histogram bins for potential %s" % P.Name)
        
        #OTHER INITIALIZATIONS
        #backtracking currently?
        self.Backtracking = False
        #first iteration
        self.StartIter = 0
        #flag for last search was successful line search (if so, then a conjugate
        #gradient direction can be picked next)
        self.PrevLineSearch = False
        
        #steps in the MD runs
        self.StepsEquil = 1000
        self.StepsProd = 100000
        self.StepsStride = 100


    def FindTarArgRanges(self, ArgDataFileName = None):
        """Finds the ranges of target argument distributions based on mins and maxs sampled."""
        #reset to original system parameters
        self.RevertSys()        
        #choose filename
        if ArgDataFileName is None:
            fn = self.ArgDataFn
        else:
            fn = ArgDataFileName
        if self.SaveLoadArgData and os.path.isfile(fn):
            print "Loading pre-compiled target argument data from %s" % fn
            fo = file(fn, 'r')
            self.RefPos, self.SavedTarArgData = pickle.load(fo)
            fo.close()
            self.ModSys.ForceField.Arg.SetData(self.SavedTarArgData) 
        else:
            base.ParseArgRanges(self.TarTraj, self.ModSys)
     

    def ParseTarData(self, AutoHistRanges = True, ArgDataFileName = None):
        """Reads target trajectory, parses argument data, and sets argument mins & maxes. 
AutoHistRanges:  True to set ranges for histograms based on argument ranges        
"""
        #reset to original system parameters
        self.RevertSys()
        #choose filename        
        if ArgDataFileName is None:
            fn = self.ArgDataFn
        else:
            fn = ArgDataFileName
                
        #parse ff component argument min, max, avg, and distributions from trajectory
        if self.SaveLoadArgData and os.path.isfile(fn):
            print "Loading pre-compiled target argument data from %s" % fn
            print "WARNING: if system settings have changed, must delete old arg data and restart."
            fo = file(fn, 'r')
            self.RefPos, self.SavedTarArgData = pickle.load(fo)
            fo.close()
            self.ModSys.ForceField.Arg.SetData(self.SavedTarArgData)  
               
        else:         
            #set the target argument ranges if they have not yet been calculated        
            if AutoHistRanges:
                self.FindTarArgRanges()
                for P in self.ModSys.ForceField:
                    P.Arg.SetupHist()
                    
            if self.Verbose: print "Initial parsing of trajectory."            
            self.RefPos, TrajRefBoxL = base.ParseArgHist(self.TarTraj, self.ModSys, self.Beta)
            if not TrajRefBoxL is None:
                if self.UseTrajBoxL:
                    self.ModSys.ScaleBox(TrajRefBoxL)
                    self.RefBoxL = TrajRefBoxL
            #reference positions
            if not self.RefPosInd is None:
                self.RefPos = self.TarTraj[self.RefPosInd]
            #save distributions
            self.SavedTarArgData = self.ModSys.ForceField.Arg.GetData()
            #dump out data for quick loading next time
            if not fn is None:
                fo = file(fn, 'w')
                pickle.dump((self.RefPos, self.SavedTarArgData), fo)
                fo.close()
            
        #set bounds and constraints for ill-sampled parameters
        if self.Verbose:
            print "Identifying bounds and constraints."
        for P in self.ModSys.ForceField:
            P.SetBounds(MaxFracChange = self.MaxFracChange)
            
        self.TarIsParsed = True
    
    
    def CheckReady(self):
        if not self.TarIsParsed:
            print "ParseTarData() not previously called.  Now calling..."
            self.ParseTarData()
        
    def Reset(self):
        """Resets optimizer for completely new iterations."""
        optimize.OptimizeClass.Reset(self)
        self.SrelShift = 0.
        self.SrelReweightShift = None
        
    def Cleanup(self):
        """Cleans up this class for deletion."""
        self.DelTarTraj()
        self.DelModTraj()

        
    def __del__(self):
        """Garbage collecting routines."""
        self.Cleanup()
        optimize.OptimizeClass.__del__(self)
        

    def __setattr__(self, name, val):
        optimize.OptimizeClass.__setattr__(self, name, val)
        if name == "FilePrefix":
            if val is None:
                self.LogBuffer = ""
                self.SumFn = None
                self.TarHistFn = None
                self.ModHistFn = None
                self.ArgDataFn = None
                self.PotentialsFn = None
            else:
                self.LogBuffer = ""
                self.LogFn = val + "_log.txt"
                self.SumFn = val + "_ff.dat"
                self.PotentialsFn = val + "_potentials.txt"                
                self.TarHistFn = val + "_hist.txt"   
                self.ModHistFn = val + "_modhist.txt"  
                self.ArgDataFn = val + "_argdata.pickle"
                if os.path.exists(self.LogFn): os.remove(self.LogFn)               
    
    
    def RevertSys(self):
        self.ModSys.TempSet = self.TempSet
        self.ModSys.ScaleBox(self.RefBoxL)
        if not self.SavedTarArgData is None:
            self.ModSys.ForceField.Arg.SetData(self.SavedTarArgData)
        if not self.RefPos is None:
            self.ModSys.Pos[:] = self.RefPos

            
    def GetSplineKnotHists(self):
        ret = []
        for P in self.ModSys.ForceField:
            if P.IsSpline:
                ret.append(P.KnotHist())
            else:
                ret.append(None)
        return ret
    
    def GetFLSplineKnotHists(self):
        ret = []
        for P in self.ModSys.ForceField:
            if P.UsesFluct:
                ret.append(P.FLKnotHist())
            else:
                ret.append(None)
        return ret
        

    def UpdateTempFilePrefix(self):
        if self.TempFilePrefix is None:
            self.TempFilePrefix = utility.NewTempFilePrefix(Dir = self.TempFileDir)
            if self.Verbose:
                print "Temporary file prefix is " + self.TempFilePrefix

    
    def UpdateModTraj(self):
        """Makes a trajectory for the model system."""  
        import sim.alg.bennett as bennett    
        #rotate old to new        
        self.RotateModTraj()
        #if we are not reweighting just get rid of the old
        if not self.ReweightOldModTraj:
            self.DelOldModTraj()
        #get a temporary file prefix
        self.UpdateTempFilePrefix()
        NewModPrefix = self.TempFilePrefix + "1"
        if not self.OldModTraj is None:
            if self.OldModTraj.FileName.startswith(NewModPrefix):
                NewModPrefix = self.TempFilePrefix + "2"
        if self.Verbose:
            print "Making model trajectory for reweighting" 
        #setup system
        self.RevertSys()   
        #make the trajectory
        self.NewModTraj = self.MakeModTrajFunc(self.ModSys, NewModPrefix, 
                                               self.StepsEquil, self.StepsProd, self.ElecSys, self.StepsStride,
                                               Verbose = self.Verbose)
        #add the parameter set so that we remember it
        self.NewModTraj.Param = self.Param.copy()
        
        if not self.ReweightOldModTraj or self.OldModTraj is None:
            #just point directly to the new trajectory
            self.ModTraj = self.NewModTraj  
        else:
            #get energies of old trajectory in new forcefield
            if self.Verbose: print "Getting cross energies of old trajectory in new forcefield."
            EneOO, EneON, LogWeightO, Nullx = base.GetEneTraj(self.OldModTraj, self.ModSys)
            #change parameters and get opposite energies
            if self.Verbose: print "Getting cross energies of new trajectory in old forcefield."
            self.SetParam(self.OldModTraj.Param)
            EneNN, EneNO, LogWeightN, Nullx = base.GetEneTraj(self.NewModTraj, self.ModSys)
            #change parameters back
            self.SetParam(self.NewModTraj.Param)
                        
            if self.Verbose: print "Running Bennet algorithm to reweight old to new trajectory."
            try:
                FE, FEerr, LogWeightN, LogWeightO = bennett.BennettFE(EneNN, EneON, EneOO, EneNO, 
                                                                      Beta1 = self.Beta, Beta2 = self.Beta,
                                                                      Reweight = True, Verbose = True)
            except bennett.BennettException:
                #just use the new trajectory
                if self.Verbose: print "Could not reweight... just using new trajectory instead."
                self.ModTraj = self.NewModTraj
            else:
                #make a combined trajectory
                if self.Verbose:
                    FracN = np.sum(np.exp(LogWeightN))
                    FracO = np.sum(np.exp(LogWeightO)) 
                    FracN, FracO = FracN / (FracO + FracN), FracO / (FracO + FracN)
                    print "Reweighting starting with %.1f%% new and %.1f%% old trajectory" % (FracN*100, FracO*100)
                NewFrameData = {"LogWeight" : [LogWeightN, LogWeightO],
                                "PEnergy" : [EneNN, EneON]}
                self.ModTraj =  traj.Multi([self.NewModTraj, self.OldModTraj],
                                           NewFrameData = NewFrameData)
                    
        #update the min reweighting fraction
        if not self.MinReweightFrames is None:
            if self.MinReweightFrames > len(self.ModTraj):
                raise ValueError("MinReweightFrames is bigger than trajectory length.")
            self.MinReweightFrac = float(self.MinReweightFrames) / len(self.ModTraj)
            
        #update the srel shifts
        self.SrelShift = self.Srel
        self.SrelReweightShift = None
        #reweighting fraction is 100% (no reweighting yet)
        self.ReweightFrac = 1.


    def DelOldModTraj(self):
        """Deletes the old model trajectory."""
        if not self.OldModTraj is None:
            print "> Deleting temporary trajectory file %s" % self.OldModTraj.FileName
            self.OldModTraj.Delete()
            self.OldModTraj = None               

     
    def RotateModTraj(self):
        """Rotates the new to the old model trajectory."""
        if not self.ModTraj is None:
            self.ModTraj.Close()
            self.ModTraj = None
        self.DelOldModTraj()            
        if not self.NewModTraj is None:
            self.OldModTraj = self.NewModTraj
            self.NewModTraj = None


    def DelModTraj(self):
        """Deletes a model trajectory, if any."""
        if not self.ModTraj is None:
            self.ModTraj.Close()
            self.ModTraj = None
        self.DelOldModTraj()
        if not self.NewModTraj is None:
            print "> Deleting temporary trajectory file %s" % self.NewModTraj.FileName
            self.NewModTraj.Delete()
            self.NewModTraj = None               

                
    def DelTarTraj(self):
        """Deletes a target trajectory, if any."""
        if not self.TarTraj is None:
            self.TarTraj.Close()
            self.TarTraj = None
            if os.path.isfile(self.TarTrajFn):
                os.remove(self.TarTrajFn)
                print "> Deleted temporary trajectory file %s" % self.TarTrajFn
                self.TarTrajFn = ""        
            
        
    def OutputLog(self, Flush = True, OtherVals = []):
        """Update the log file after an iteration."""
        #file output
        if not self.LogFn is None:
            Names = []
            for (i, n) in enumerate(self.Names):
                n = n.replace(" ", "")
                if i in self.Constrained:
                    Names.append(n.strip() + "*")
                else:
                    Names.append(n)
            sMode = "".join([x[0] for x in self.Mode.split()])
            s = ""
            LogExists = os.path.isfile(self.LogFn)
            if not LogExists:
                #write the header
                s += "%-10s %-4s %-12s %-12s %-12s %-12s " % ("Iter", "Mode", "Srel", "Bias", "Obj", "dx")               
                for (nm, val) in OtherVals:
                    s += "%-12s " % nm
                s += " ".join(["%-20s" % n for n in Names])
                s += "\n"                                 
            #output parameter values
            s += "%-10d %-4s %-12.4e %-12.4e %-12.4e %-12.4e " % (self.Iter, sMode, self.Srel, 
                                                                  self.Bias, self.Obj, self.dx)
            for (nm, val) in OtherVals:
                s += "%-12.4e " % val            
            s += " ".join(["%-20.4e" % p for p in self.Param])
            s += "\n"
            self.LogBuffer += s
            if Flush:
                if LogExists:
                    Mode = "a"
                else:
                    Mode = "w"
                try:
                    file(self.LogFn, Mode).write(self.LogBuffer)
                except IOError:
                    print "WARNING: Could not write to log file.  Will try again next time."
                else:
                    self.LogBuffer = ""
       


    def OutputSumFile(self, Param = None):
        """Update the summary file after an iteration.
Param: if supplied, use this parameter set rather than the current one."""
        if not Param is None:
            OldParam = self.ModSys.ForceField.Param.Val.copy()
            self.ModSys.ForceField.SetParam(Param)
        if not self.SumFn is None:
            s = self.ModSys.ForceField.ParamString()
            file(self.SumFn, "w").write(s) 
        if not self.PotentialsFn is None:
            s = self.OutputReport()
            file(self.PotentialsFn, "w").write(s) 
        if not Param is None:
            self.ModSys.ForceField.SetParam(OldParam)                

    def GetTarHists(self, HistLabel = "target"):
        """Gets a list of target histogram objects."""
        self.RevertSys()
        #get the target list of histograms
        Hists = []
        for P in self.ModSys.ForceField:
            Hists.extend(P.HistList(HistLabel = HistLabel))
        return Hists

    def OutputTarHistFile(self):
        """Update the target histogram file."""
        if not self.TarHistFn is None:
            Hists = self.GetTarHists()
            s = "\n========TARGET ARGUMENT HISTOGRAMS========\n"
            s = "\n".join([h.HistString() for h in Hists])
            file(self.TarHistFn, "w").write(s)
        
    def GetModHists(self, HistLabel = "model"):
        """Gets a list of model histogram objects from current model trajectory."""
        if self.ModTraj is None:
            raise ValueError("Must make a model trajectory first before getting histograms.")
        #revert system
        self.RevertSys()
        #parse the model trajectory into a new histogram
        base.ParseArgHist(self.ModTraj, self.ModSys, self.Beta, Reweight = True)
        #get the model histograms
        Hists = []
        for P in self.ModSys.ForceField:
            Hists.extend(P.HistList(HistLabel = HistLabel))
        #restore old (target) histogram data
        self.RevertSys()
        #return
        return Hists
        
    def GetBothHists(self, TarHistLabel = "target", ModHistLabel = "model"):
        """Returns target plus model hists."""
        Hists = []
        TarHists = self.GetTarHists(HistLabel = TarHistLabel)
        ModHists = self.GetModHists(HistLabel = ModHistLabel)
        for (tHist, mHist) in zip(TarHists, ModHists):
            Hists.append(tHist + mHist)
        return Hists
            
    def OutputModHistFile(self):
        """Updates the model histogram file."""
        if not self.ModHistFn is None:
            Hists = self.GetBothHists()
            #output to a file
            s = "\n========TARGET AND MODEL ARGUMENT HISTOGRAMS========\n"
            s += "\n".join([h.HistString() for h in Hists])
            file(self.ModHistFn, "w").write(s)
            
    def OutputPlot(self, UseIter = False, FilePrefix = None):
        """Updates the plot of histograms and distributions."""
        if FilePrefix is None:
            FilePrefix = self.FilePrefix
        if FilePrefix is None or not EnablePlots:
            return
        import plot_forcefield
        if not self.ModHistFn is None:
            if UseIter:
                PlotFn = "%s_potentials%06d.%s" % (FilePrefix, self.Iter, PlotFmt)
            else:
                PlotFn = "%s_potentials.%s" % (FilePrefix, PlotFmt)
            plot_forcefield.PlotForceField(self.ModSys, PlotFn, self.ModHistFn)
        #output the time progression
        if not self.LogFn is None and self.Iter > self.StartIter:
            PlotFn = "%s_iterations.%s" % (FilePrefix, PlotFmt)
            plot_forcefield.PlotParameterHistory(self.LogFn, PlotFn)
            
    
    def OutputPotentials(self):
        """Outputs potentials."""
        s = ""
        for P in self.ModSys.ForceField:
            s += P.ReportString() + "\n"
        file(self.PotentialsFn, "w").write(s)
           

    def CalcObjTarTraj(self, CalcDeriv = True):
        """Computes the objective function from the trajectory entirely.
CalcDeriv: computes derivatives wrt ff parameters if True"""
        if not self.TarIsParsed:
            raise ValueError("Must run ParseTarData first.")
        #setup system
        self.RevertSys()
        #calc srel
        ret = base.CalcSrelTraj(self.TarTraj, self.ModSys, self.Beta, CalcDeriv = CalcDeriv, 
                                Penalties = self.Penalties, CalcFluct = self.UseFluct)
        self.Obj, self.DObj, self.DDObj, AvgFluctTerm = ret
        self.Srel = self.Obj
        self.Bias = 0.
        if self.UseFluct:
            self.FluctBias = AvgFluctTerm
            self.Obj += AvgFluctTerm
        #set the reweight fraction to 1 since it doesn't have much meaning here
        self.ReweightFrac = 1.
        #calculate objective and compensate for constraints
        self.UpdateObj()  
        #return the objective
        return self.Obj
            
            
    def CalcObjModTraj(self, CalcDeriv = True):
        """Computes the objective function and derivatives using a recorded trajectory from an integrator."""
        if not self.TarIsParsed:
            raise ValueError("Must run ParseTarData first.")
        #set reweighting fraction to 1
        self.ReweightFrac = 1.
        #setup system
        self.RevertSys()

        #get the target averages
        if self.UseTarHists:
            TarDUParam, TarDDUParam, TarPEnergy = base.CalcAvgsHist(self.ModSys, CalcDeriv = CalcDeriv)
            AvgFluctTerm = 0.
        else:
            ret = base.CalcAvgsTraj(self.TarTraj, self.ModSys, self.Beta, 
                                    Reweight = False, CalcDeriv = CalcDeriv, CalcFluct = self.UseFluct)
            TarDUParam, TarSqDUParam, TarDDUParam, TarCount, ReweightFrac, SrelReweightTerm, TarPEnergy, AvgFluctTerm = ret
            
        #get the model averages
        ret = base.CalcAvgsTraj(self.ModTraj, self.ModSys, self.Beta, Reweight = True, 
                                CalcDeriv = CalcDeriv, Penalties = self.Penalties)

        ModDUParam, ModSqDUParam, ModDDUParam, ModCount, ReweightFrac, SrelReweightTerm, ModPEnergy, Nullx = ret
        self.ReweightFrac = ReweightFrac

        if base.DEBUG:
            print "%20s %12s %12s %12s %12s %s" % ("Param", "B*<dU/dl>T", "B*<dU/dl>M", "B2*<dU/dl2>M", "B2*<dU/dl>2M", "ModCount")
            for (i,(t,m,m2)) in enumerate(zip(TarDUParam, ModDUParam, ModSqDUParam.diagonal())):
                if self.Fixed[i]: continue
                print "%-20s %12.4e %12.4e %12.4e %12.4e %10.2e" % (self.Names[i], 
                                                                    self.Beta * t,
                                                                    self.Beta * m, 
                                                                    self.Beta**2 * m2,
                                                                    self.Beta**2 * m*m,
                                                                    ModCount[i])
        #compute the Srel derivatives
        self.DObj = self.Beta * (TarDUParam - ModDUParam)
        self.DDObj = self.Beta * (TarDDUParam - ModDDUParam 
                                  + self.Beta * ModSqDUParam 
                                  - self.Beta * np.outer(ModDUParam, ModDUParam))
        #set the derivative for fixed parameters to zero
        self.DObj[self.Fixed] = 0.
        self.DDObj[np.ix_(self.Fixed, self.Fixed)] = 0.
    
        #compute srel
        if self.SrelReweightShift is None:
            self.SrelReweightShift = self.Beta * TarPEnergy + SrelReweightTerm
        self.Srel = self.SrelShift + self.Beta * TarPEnergy + SrelReweightTerm - self.SrelReweightShift
            
        self.Obj = self.Srel + AvgFluctTerm
        self.FluctBias = AvgFluctTerm
        self.Bias = 0.
                
        #calculate objective and compensate for constraints 
        self.UpdateObj() 
                   
        #return the objective
        return self.Obj
            
    
    def CalcObj(self, CalcDeriv = True):
        """Computes the objective function and derivatives using current class options."""
        if self.ModTraj is None:
            raise AttributeError("No model trajectory has been made yet.")        
        if self.ReweightTar:
            return self.CalcObjTarTraj(CalcDeriv)
        else:
            return self.CalcObjModTraj(CalcDeriv)
        

    def TestObj(self, SkipZeroDObj = True, DoDDObj = False, 
                DoDUParam = False, DoDWParam = False,
                ParamDelta = 0.00001, ParamMult = 1.000001):
        """Computes Obj derivatives from finite difference and compares
to explicitly calculated."""
        self.CheckReady()            
        print "Test of DObj calculation:"
        #initialize any constraints
        self.InitConstraints(Verbose = self.Verbose)
        self.SetParam(self.Param)
        #define shortcut objective and derivative getting function
        def GetObj():
            self.CalcObj()
            return self.Obj, self.DObj.copy(), self.DDObj.copy()
        #setup system
        self.RevertSys()
        #finite differences
        Obj1, DObj1, DDObj1 = GetObj()
        for (i, Pen) in enumerate(self.Penalties):
            print "Penalty #%d %s" % (i, Pen.Name)
            print "   average: %12.5e" % Pen.Avg
            print "   target:  %12.5e" % Pen.Target
            print "   bias:    %12.5e" % Pen.Obj
        Param1 = self.Param.copy()
        for (i, iFixed) in enumerate(self.Fixed):
            if iFixed or i in self.Constrained or (SkipZeroDObj and DObj1[i]==0): 
                continue
            iVal1 = Param1[i]
            iDObj1 = DObj1[i]
            if iVal1 == 0.:
                iVal2 = ParamDelta
            else:
                iVal2 = iVal1 * ParamMult
            Param2 = Param1.copy()
            Param2[i] = iVal2
            self.SetParam(Param2)
            Obj2, DObj2, DDObj2 = GetObj()
            iDObj2 = DObj2[i]
            iVal0 = 2*iVal1 - iVal2
            Param0 = Param1.copy()
            Param0[i] = iVal0
            self.SetParam(Param0)
            Obj0, DObj0, DDObj0 = GetObj()
            iDObj0 = DObj0[i]
            DObjDiff = (Obj2 - Obj0) / (iVal2 - iVal0)
            print "Param %s:" % self.Names[i]
            print "   DObjComp: %12.5e" % iDObj1
            print "   DObjDiff: %12.5e  err: %12.5e" % (DObjDiff, abs(iDObj1 - DObjDiff))
            if not DoDDObj: continue
            for (j, jFixed) in enumerate(self.Fixed):
                if jFixed or j in self.Constrained or (SkipZeroDObj and DObj1[j]==0):
                    continue
                if i == j:
                    DDObjDiff = (Obj2 - 2*Obj1 + Obj0) / ((iVal2 - iVal1)*(iVal1- iVal0))
                    DDObjDiff2 = (iDObj2 - iDObj0) / (iVal2 - iVal0)
                    DDObjComp = DDObj1[i,j]
                else:
                    ijDDObj1 = DDObj1[i,j]
                    ijDDObj2 = DDObj2[i,j]
                    jVal1 = Param1[j]
                    if jVal1 == 0.:
                        jVal2 = ParamDelta
                    else:
                        jVal2 = jVal1 * ParamMult
                    Param3 = Param1.copy()
                    Param3[j] = jVal2
                    self.SetParam(Param3)
                    Obj3 , DObj3, DDObj3 = GetObj()
                    iDObj3 = DObj3[i]
                    ijDDObj3 = DDObj3[i,j]
                    Param4 = Param2.copy()
                    Param4[j] = jVal2
                    self.SetParam(Param4)
                    Obj4, DObj4, DDObj4 = GetObj()
                    iDObj4 = DObj4[i]
                    ijDDObj4 = DDObj4[i,j]
                    DDObjComp = np.mean([ijDDObj1, ijDDObj2, ijDDObj3, ijDDObj4])
                    DDObjDiff = (Obj4 - Obj3 - Obj2 + Obj1)/(iVal2 - iVal1)/(jVal2 - jVal1)
                    DDObjDiff2 = (iDObj3 + iDObj4 - iDObj1 - iDObj2) / (jVal2 - jVal1) / 2.
                print "Params %s and %s:" % (self.Names[i], self.Names[j])
                print "  DDObjComp: %12.5e" % DDObjComp
                #print "  DDObjDiff: %12.5e  err: %12.5e" % (DDObjDiff, abs(DDObjComp - DDObjDiff))
                print "  DDObjDiff: %12.5e  err: %12.5e" % (DDObjDiff2, abs(DDObjComp - DDObjDiff2))
        self.SetParam(Param1)   
        
        if DoDUParam: 
            print "-"*40        
            self.ModSys.Flags.CalcDUParam = True
            self.ModSys.ForceField.Eval()
            DUParam1 = self.ModSys.ForceField.Param.DU
            for (i, iFixed) in enumerate(self.Fixed):
                if iFixed or i in self.Constrained: 
                    continue
                iVal1 = Param1[i]
                iDUParam1 = DUParam1[i]
                if iVal1 == 0.:
                    iVal2 = ParamDelta
                else:
                    iVal2 = iVal1 * ParamMult
                Param2 = Param1.copy()
                Param2[i] = iVal2
                self.SetParam(Param2)
                self.ModSys.ForceField.Eval()            
                U2 = self.ModSys.PEnergy
                iVal0 = 2*iVal1 - iVal2
                Param0 = Param1.copy()
                Param0[i] = iVal0
                self.SetParam(Param0)
                self.ModSys.ForceField.Eval()
                U0 = self.ModSys.PEnergy
                DUParamDiff = (U2 - U0) / (iVal2 - iVal0)
                print "Param %s:" % self.Names[i]
                print "   DUParamComp: %12.5e" % iDUParam1
                print "   DUParamDiff: %12.5e  err: %12.5e" % (DUParamDiff, abs(iDUParam1 - DUParamDiff)) 
                self.SetParam(Param1)    
                
        if DoDWParam:
            print "*"*40
            self.ModSys.Flags.CalcVirial = True 
            self.ModSys.Flags.CalcDWParam = True 
            self.ModSys.ForceField.Eval()
            DWParam1 = self.ModSys.ForceField.Param.DW 
            for (i, iFixed) in enumerate(self.Fixed):
                if iFixed or i in self.Constrained: 
                    continue
                iVal1 = Param1[i]
                iDWParam1 = DWParam1[i]
                if iVal1 == 0.:
                    iVal2 = ParamDelta
                else:
                    iVal2 = iVal1 * ParamMult
                Param2 = Param1.copy()
                Param2[i] = iVal2
                self.SetParam(Param2)
                self.ModSys.ForceField.Eval()            
                W2 = self.ModSys.Virial 
                iVal0 = 2*iVal1 - iVal2
                Param0 = Param1.copy()
                Param0[i] = iVal0
                self.SetParam(Param0)
                self.ModSys.ForceField.Eval()
                W0 = self.ModSys.Virial 
                DWParamDiff = (W2 - W0) / (iVal2 - iVal0)
                print "Param %s:" % self.Names[i]
                print "   DWParamComp: %12.5e" % iDWParam1
                print "   DWParamDiff: %12.5e  err: %12.5e" % (DWParamDiff, abs(iDWParam1 - DWParamDiff)) 
                self.SetParam(Param1)             
        
           
    def TestObjAll(self, *args, **kwargs):
        """Runs tests separately on penalties."""
        self.UpdateMode = 1
        print "======== REPORT: derivatives computed for Srel ========"
        self.TestObj(*args, **kwargs)
        self.UpdateMode = 2
        print "======== REPORT: derivatives computed for Bias ========"
        self.TestObj(*args, **kwargs)
        Coefs = [Pen.Coef for Pen in self.Penalties]
        for (iPen, iCoef) in zip(self.Penalties, Coefs):
            for jPen in self.Penalties:
                jPen.Coef = 0
            iPen.Coef = iCoef
            print "======== REPORT: derivatives computed for penalty %s ========" % iPen.Name
            self.TestObj(*args, **kwargs)
            for (jPen, jCoef) in zip(self.Penalties, Coefs):
                jPen.Coef = jCoef
        self.UpdateMode = 0
            
                    
    def DirModeStr(self):
        """Returns a string describing the current search direction mode."""
        s = DirModes[self.DirMode]
        if self.Backtracking:
            s += " -BACKTRACK"
        return s            

    
    def Initdx(self):
        """Resets dx values to defaults."""
        self.dxMode = dxModeDefaults.copy()
        
    def Setdx(self, ScaleFactor = None):
        """Sets and saves the current value of dx."""
        #scale if requested
        if not ScaleFactor is None:
            self.dx = self.dx * ScaleFactor
        #clip to max bounds
        self.dx = min(self.dx, dxMaxMode[self.DirMode])
        #save the dx value for this mode
        self.dxMode[self.DirMode] = self.dx
        
    def Getdx(self):
        """Gets an initial estimate of dx for this mode."""
        self.dx = self.dxMode[self.DirMode]
        

    def CalcSearchDir(self, ForceSteepest = False):
        """Returns a search direction."""
          
        #record last direction
        self.LastDir = self.Dir.copy()
        
        #find which parameters to use for computing direction
        Ignore = self.Fixed.copy()
        for i in self.Constrained:
            Ignore[i] = True
        #ignore any parameters with zero hessian
        for (i, ThisDObj) in enumerate(self.DDObj):
            if all(self.DDObj[i,:] == 0): 
                Ignore[i] = True
        Use = np.logical_not(Ignore)

        #filter parameters & gradients
        DObj = self.DObj[Use]
        DDObj = self.DDObj[np.ix_(Use, Use)] + self.HessianPad * np.identity(len(DObj))
                                
        #STEEPEST DESCENT DIRECTION
        SDDir = -DObj
        dSDDir = np.dot(SDDir, DObj)
        
        #CONJUGATE GRADIENT DIRECTION
        CGDir = None
        dCGDir = 1.e300
        self.CGBeta = 0.
        if self.PrevLineSearch and not self.LastDObj is None:
            LastDObj = self.LastDObj[Use]
            LastDir = self.LastDir[Use]
            #choose an appropriate way to determine Beta
            self.CGRestartCount += 1
            if not self.DirMode in ["SD", "CG"]:
                #Beta should be zero because previous step was not a CG step
                b = 0
            elif self.CGRestartCount >= self.CGRestartIters:
                #restart the search with Beta=0
                self.CGRestartCount = 0
                b = 0
            elif self.CGFormula == 0:
                #Fletcher-Reeves
                b = np.sum(DObj**2) / np.sum(LastDObj**2)
            elif self.CGFormula == 1:
                #Polak-Ribiere
                b = np.sum(DObj * (DObj - LastDObj)) / np.sum(LastDObj**2)
                b = max(0., b)
            else:
                raise ValueError("Did not recognize CGFormula == %d" % self.CGFormula)
            CGDir = -DObj + b * LastDir
            dCGDir = np.dot(CGDir, DObj)
            self.CGBeta = b
            
        #HESSIAN EIGEN DIRECTION
        #get the hessian direction, from Newton-Raphson equation H^-1 * Grad
        HDir = np.zeros_like(DObj)
        dHDir = 1.e300
        if self.UseHessian and len(DObj) > 1:
            #eigendecompose
            try:
                eigs, vecs = np.linalg.eigh(DDObj)
            except np.linalg.LinAlgError:
                eigs = [-1]
                vecs = None
            #if there are any negative eigs we shouldn't use this
            if np.all(eigs >= 0):
                #rebuild the hessian inverse, optionally filtering eigenvalues
                dHDir = 0.
                Delta = np.zeros_like(HDir)
                for (i, eig) in enumerate(eigs):
                    if eig < self.MinEig: continue
                    vec = vecs[:,i]
                    coef = np.dot(vec, DObj)
                    #find the direction contribution from this eigenvector
                    d = -coef * vec / eig
                    #find the expected change in the objective function
                    Delta[i] = -coef*coef / eig
                    #include the direction 
                    #check if the magnitude of the direction is too big
                    if np.max(np.abs(d)) <= self.HessianMaxDir:
                        HDir = HDir + d
                        dHDir = dHDir + Delta[i]
                    elif base.DEBUG or DEBUGHESSIAN:
                        print "Omitting direction:\n", d
            #save the hessian
            if (base.DEBUG or DEBUGHESSIAN) and not vecs is None:
                ThisNames = [self.Names[i] for (i,u) in enumerate(Use) if u]
                sThisNames = " ".join(ThisNames)
                n = len(DObj)
                Contrib = [-np.dot(vecs[:,i], DObj) / eigs[i] for i in range(n)]
                Contrib = np.array(Contrib) 
                rContrib = np.abs(Contrib) / np.sum(np.abs(Contrib))
                s = "HESSIAN\n"
                s += sThisNames + "\n"
                for i in range(n):
                    s += " ".join(["%15.4e" % DDObj[i,j] for j in range(n)]) + "\n"
                s += "\nEIGENANALYIS\n"
                s += "%15s " % ("eig",) + sThisNames + "\n"
                for i in range(n):
                    s += "%15.4e " % (eigs[i],)
                    s += " ".join(["%15.4e" % vecs[j,i] for j in range(n)]) + "\n"
                s += "\n%15s %15s %15s %15s\n" % ("eig", "contrib_to_dir", "rel_contrib", "max_param_contrib")
                for i in range(n):
                    ind = np.argmax(np.abs(vecs[:,i]))
                    s += "%15.4e %15.4e %8.2f%% %15s\n" % (eigs[i], Contrib[i], rContrib[i]*100, ThisNames[ind])                  
                s += "\nHDIR\n"
                s += "%15s %15s %15s\n" % ("Param", "DObj", "HDir")
                for i in range(n):
                    s += "%15s %15.4e %15.4e\n" % (ThisNames[i], DObj[i], HDir[i])
                file(self.FilePrefix + "_hessian.txt", "w").write(s)                

        #CHOOSE WHICH DIRECTION TO USE
        self.Dir.fill(0.)
        if ForceSteepest or (dHDir > 0 and dCGDir > 0):
            #use steepest descent if hessian or CG directions take us uphill
            #or if we have taken too many CG steps in a row
            self.Dir[Use] = SDDir
            self.DirMode = "SD"
        elif (dHDir < 0 and self.UseHessian):
            #use hessian descent if it takes us downhill or if we
            #are not to use conjugate gradient
            self.Dir[Use] = HDir
            self.DirMode = "HD"
        else:
            #use conjugate gradient direction if otherwise
            self.Dir[Use] = CGDir
            self.DirMode = "CG"
            

    def LineSearch(self):
        """Performs a line search to minimize Srel.
Returns False if not successful such that reweighting needs to happen."""       
        #get the starting parameters
        StartParam = self.Param.copy()
        ThisIter = 0
        MaxIter = self.LSMaxIter

        self.Mode = "LINE SEARCH BRACKET"
        self.PrevLineSearch = False
        self.Backtracking = False
               
        #make shortcut functions for evaluating the objective
        EvalList = [] 
        def SetParamDist(dist):
            NewParam = StartParam + dist * self.Dir
            self.SetParam(NewParam)   
        
        def EvalObj(dist):
            self.LastParam = self.Param.copy()
            SetParamDist(dist)
            self.CalcObj(CalcDeriv = False)
            self.OutputLog(Flush = False)
            if self.Verbose:
                print self.Output1(ShowParams = self.LSVerbose)
            EvalList.append((dist, self.Obj)) 
            return self.Obj
            
        def ReweightFlag():
            #returns true if reweighted trajectory fraction is too low
            return np.any(self.ReweightFrac < self.MinReweightFrac)
            
            
        #find which variable will exceed range first
        DistLim, SoftDistLim = self.CalcDistLim()
        
        #get an estimate for dx
        self.Getdx()

        #make sure dx is small enough not to cause go past boundary
        ct = 0
        while self.dx >= DistLim:
            self.dx = self.dx * self.LSBackScale
            ct += 1
            if ct > 10000:
                raise ValueError("Error in backscaling: DistLim = %.4e, dx = %.4e" % (DistLim, self.dx))

        #take the first two steps and compute energies
        Obj0 = self.ObjVals[-1]
        ReweightFracs = [np.min(self.ReweightFrac)]
        Dists = [0., self.dx]
        Vals = [Obj0, EvalObj(Dists[1])]
        ReweightFracs += [np.min(self.ReweightFrac)]

        #if the second point is the same, try increasing step size
        while Vals[0] == Vals[1]:
            if self.Verbose:
                print "Step size looks too small; increasing\n"
            if Dists[1] >= DistLim:
                print "ERROR: No change in objective function with dx = %11.4e and exceeding distance limit." % self.dx
                raise ValueError("No change in objective function.")
                return
            else:
                self.dx = self.dx * self.LSAccel
            Dists[1] = self.dx
            Vals[1] = EvalObj(self.dx)
            ReweightFracs[1] = np.min(self.ReweightFrac)
            ThisIter += 1
            if ThisIter > MaxIter:
                print "ERROR: Could not produce a change in objective function with dx = %11.4e" % self.dx
                raise ValueError("No change in objective function.")
                return
     
        #if the second point is not downhill, back
        #off and take a shorter step until we find one;
        #also back off if reweight fraction is too small
        while True:        
            if np.isnan(Vals[1]):
                if self.Verbose:
                    print "Step size too large: got NaN.  Backing off.\n"
            elif Vals[1] > Vals[0]:
                if self.Verbose:
                    print "Step size too large: objective increased.  Backing off.\n"
            elif ReweightFlag():
                if self.Verbose:
                    print "Step size too large: surpassed reweighting limit.  Backing off.\n"
            elif ReweightFracs[0] - ReweightFracs[1] > self.MaxDeltaReweightFrac:
                if self.Verbose:
                    print "Step size too large: reweight fraction decreased too fast.  Backing off.\n"
            else:
                break
            self.dx = self.dx * self.LSBackScale
            self.Backtracking = True
            Dists[1] = self.dx
            Vals[1] = EvalObj(self.dx)
            ReweightFracs[1] = np.min(self.ReweightFrac)
            ThisIter += 1
            if ThisIter > MaxIter:
                print "ERROR: could not find a step size that reduced objective function."
                raise ValueError("Could not decrease objective function.")
                return                
                   
        #save the current estimate for a good dx
        self.Setdx()

        if self.Verbose:
            print "Finding three points to bracket minimum\n"  
               
        #keep stepping forward until the third point is higher;
        #then we have bracketed a minimum
        while Vals[-1] <= Vals[-2]:
            ThisIter += 1
            self.Backtracking = False
            #check limits
            if self.dx + Dists[-1] >= DistLim:
                self.dx = (DistLim - Dists[-1]) * 0.5
            #find a next point and evaluate
            Dists = Dists + [Dists[-1] + self.dx]               
            Vals = Vals + [EvalObj(Dists[-1])]
            ReweightFracs = ReweightFracs + [np.min(self.ReweightFrac)]
            if ReweightFlag(): 
                self.Backtracking = True
                while ReweightFlag():                  
                    #linearly interpolate to find where we think we would have reached
                    #the reweighting limit
                    Dist0 = Dists[-2] + (self.MinReweightFrac - ReweightFracs[-2]) * \
                            (Dists[-1] - Dists[-2]) / (ReweightFracs[-1] - ReweightFracs[-2])
                    #now go 2/3 of the way there
                    self.dx = Dists[-1] - Dists[-2]
                    self.Setdx()    
                    Dists[-1] = 0.33 * Dists[-2] + 0.67 * Dist0    
                    Vals[-1] = EvalObj(Dists[-1])
                self.Backtracking = False
                return False
            #check iterations
            if ThisIter > MaxIter:
                print "Exceeded maximum number of line search steps (%d)." % MaxIter
                return True
                 
            if len(Vals) == 3: continue        
            
            #check if we increased too much; if so, back off
            if Vals[-1] - Vals[-2] > self.LSMaxIncr * abs(Vals[-2]):
                Vals = Vals[:3]
                Dists = Dists[:3]
                ReweightFracs = ReweightFracs[:3]
                self.dx = self.dx * self.LSBackScale
                self.Backtracking = True
            else:
                #shift all of the points over
                Vals = Vals[-3:]
                Dists = Dists[-3:]
                ReweightFracs = ReweightFracs[-3:]
                #accelerate or slow as needed
                DeltaReweightFrac = ReweightFracs[-2] - ReweightFracs[-1]
                if DeltaReweightFrac > self.MaxDeltaReweightFrac:
                    self.dx = self.dx * self.LSSlow
                elif DeltaReweightFrac < self.MinDeltaReweightFrac:
                    self.dx = self.dx * self.LSAccel
                #let's assume this is a reasonably good value of dx
                self.Setdx()

        #update an estimate for a good value of dx
        self.Backtracking = False
        self.dx = Dists[1] - Dists[0]        
        self.Setdx()

        if self.Verbose:
            print "Bracketed minimum\n"        
                
        #we've bracketed a minimum; now we want to find it to high accuracy
        OldVal3 = Vals[1]
        UsePara = False
        ConsecPara = 0
        while True:      
            #update and check number of tries
            ThisIter += 1
            if ThisIter > MaxIter:  
                print "ERROR: Could not converge to minimum with parabolic / golden search within max iters."
                print "Showing iterations, distance, objective function."
                for (i, (dist, Obj)) in enumerate(EvalList):
                    print i, dist, Obj
                print "Dot product of DObj and Dir: ", np.sum(self.Dir, self.DObj)
                return True
            
            #store distances for ease of code-reading
            d0, d1, d2 = Dists
            Val0, Val1, Val2 = Vals
            d10 = d0 - d1
            d12 = d2 - d1

            #use a parobolic approximation to estimate the location
            #of the minimum
            Num = d12*d12*(Val0-Val1) - d10*d10*(Val2-Val1)
            Den = d12*(Val0-Val1) - d10*(Val2-Val1)
            if Den == 0:
                #parabolic extrapolation won't work; set new dist = 0 
                d3Para = 0
            else:
                #location of parabolic minimum
                d3Para = d1 + 0.5 * Num / Den

            #compute the new distance by golden section search
            Golden = 0.381966
            if d12 > d10:
                d3Gold = d1 + Golden * d12
            else:
                d3Gold = d1 - Golden * d10

            #if the parabolic estimate is out of bounds or we
            #used parabolic last two times, use golden instead
            UsePara = (d3Para < d2 and d3Para > d0) and (ConsecPara < 2)
            if UsePara:
                ConsecPara += 1
                d3 = d3Para
                self.Mode = "LINE SEARCH PARABOLIC"
            else:
                ConsecPara = 0
                d3 = d3Gold
                self.Mode = "LINE SEARCH GOLDEN"
                
            #compute the new objective
            Val3 = EvalObj(d3)            
                
            #decide which three points to keep; we want to keep
            #the three that are closest to the minimum
            if d3 < d1:
                if Val3 < Val1:
                    #get rid of point 2
                    Dists, Vals = [d0, d3, d1], [Val0, Val3, Val1]
                else:
                    #get rid of point 0
                    Dists, Vals = [d3, d1, d2], [Val3, Val1, Val2]
            else:
                if Val3 < Val1:
                    #get rid of point 0
                    Dists, Vals = [d1, d3, d2], [Val1, Val3, Val2]
                else:
                    #get rid of point 2
                    Dists, Vals = [d0, d1, d3], [Val0, Val1, Val3]
                    
            #check how much we've changed
            Val3Change = abs(OldVal3 - Val3)
            ValsDelta = abs(Vals[2] - Vals[0])
            if Val3Change < self.LSFracTol * Val3 \
               or ValsDelta < self.LSFracTol * Val3 \
               or Val3Change < self.LSAbsTol \
               or ValsDelta < self.LSAbsTol \
               or (Val3 < self.LSAbsTol and self.ReweightTar):
                #the fractional change is less than the tolerance,
                #so we are done and can exit the loop
                break
            OldVal3 = Val3

        if self.Verbose:
            print "Finished line search\n"    

        #return the at the minimum (point 1)
        SetParamDist(Dists[1])
        self.Obj = Val1
        self.PrevLineSearch = True
        
        return True
    
    
    def RunConjugateGradient(self, MaxIter = None, MaxTime = None, SteepestIter = 0, 
                             CalcHists = 2, HistsIter = 2):
        """Optimizes the relative entropy using a conjugate gradient algorithm
with trajectory reweighting.
MaxIter: number of iterations to do, if specified
MaxTime: maximum time in seconds to perform optimization
SteepestIter: force this number of first search directions to be steepest descent
CalcHists: at the end of the run how to output histograms: 0 to do nothing, 
           1 to compute histograms anew, 2 from reweighing
HistsIter: save histograms and potential figures every this iterations, 
           use a negative number to append iter to plot filename       
           """
        #set initial time and parameters
        self.CheckReady()
        StartTime = time.time()
        self.StartIter = self.Iter
        SearchDirIter = 0
        self.Mode = "INIT"
        self.Backtracking = False
        self.UseMaxChange = False
        self.PrevLineSearch = False
        self.CGRestartCount = 0
        self.Initdx()
        
        if self.Verbose:
            print "\n======== STARTING NEW CONJUGATE GRADIENT MINIMIZATION ========\n" 
            
        #initialize any constraints
        self.InitConstraints(Verbose = self.Verbose)      
        self.SetParam(self.Param)     
        
        if self.Verbose:
            print self.Output0()
        self.OutputTarHistFile()
        
        #make model trajectory
        if not self.ReweightTar:
            self.UpdateModTraj()
            self.OutputModHistFile()
            self.OutputPlot()

        #store this mask
        NotFixed = np.logical_not(self.Fixed)

        #setup for minimization
        if self.Verbose:
            print "Starting minimization\n"

        #make a shortcut function for evaluating the objective
        def EvalObj():
            self.CalcObj(CalcDeriv = True)
            if self.Verbose:
                print "Evaluating new gradient\n"

        #start conjugate gradient search
        EvalObj()
        self.CalcSearchDir(ForceSteepest = SearchDirIter < SteepestIter)
        self.Mode = self.DirModeStr()
        self.OutputLog()
        if self.Verbose:
            print self.Output2()

        #run convergence loop        
        Converged = False
        LastHistsIter = self.StartIter
        while not Converged and (MaxIter is None or MaxIter > self.Iter - self.StartIter) \
              and (MaxTime is None or time.time() - StartTime < MaxTime): 

            self.PrepUpdate()
            #do a line search
            MakeNewTrajectory = not self.LineSearch()
            
            #check if we need to make a new trajectory
            if MakeNewTrajectory and not self.ReweightTar:

                #moved too far away from the saved trajectory, make a new one
                print "Reweighting went too far; making new trajectory."  % np.min(self.ReweightFrac)
                self.UpdateModTraj()
                #evaluate the new derivatives
                EvalObj()
                Converged = False

            else:

                #evaluate the new derivatives
                EvalObj()
                #check the derivative and the tolerances
                FracTol = abs(self.Obj - self.LastObj) / abs(self.Obj)
                Test1 = FracTol < self.CGFracTol or (self.Obj < self.CGAbsTol and self.ReweightTar)
                Test2 = np.all(np.abs((self.DObj*self.Scale)[NotFixed]) < self.CGGradTol)
                if Test1:
                    print('Pass Test1: FracTol < self.CGFracTol or (self.Obj < self.CGAbsTol and self.ReweightTar)')
                if Test2:
                    print('Pass Test2: np.all(np.abs((self.DObj*self.Scale)[NotFixed]) < self.CGGradTol)')
                if (Test1 and Test2) or np.sum(NotFixed) == 1:
                    Converged = True

            #new search direction
            SearchDirIter += 1
            self.CalcSearchDir(ForceSteepest = SearchDirIter < SteepestIter)
            self.Mode = self.DirModeStr()
            self.OutputLog()
            self.OutputSumFile()
            if self.Verbose:
                print self.Output2(OtherData = [("CGBeta", self.CGBeta)])
                    
            #output detailed histograms and potentials this iteration?
            if CalcHists > 0 and self.Iter - LastHistsIter >= abs(HistsIter):
                LastHistsIter = self.Iter
                self.OutputModHistFile()
                self.OutputPlot(UseIter = HistsIter < 0)
           
        if self.Verbose:
            Elapsed = (time.time() - StartTime) / 60.
            print "Elapsed time: %.1f min" % Elapsed
            if not Converged:
                print "Minimization did not converge in %d iters\n" % self.Iter
            else:
                print "Minimization converged.\n"

        self.Mode = "FINAL"
        self.OutputLog()
        self.OutputSumFile()
        if self.Verbose:
            print self.Output3()    
            
        #output the final histograms
        if CalcHists == 1:
            self.UpdateModTraj()
            self.OutputModHistFile()
            self.OutputPlot(UseIter = HistsIter < 0)
        elif CalcHists == 2 and not self.ReweightTar:
            self.OutputModHistFile()
            self.OutputPlot(UseIter = HistsIter < 0)

        return Converged            


    def RunNewtonRaphson(self, MinIter = None, MaxIter = None, SteepestIter = 0,
                         MaxTime = None, CalcHists = 1, HistsIter = 20):
        """Optimizes the relative entropy using an iterative Newton Raphson technique.
MinIter: minimum number of iterations to do, if specified
MaxIter: maximum number of iterations to do, if specified
SteepestIter: force this number of first iterations to be steepest descent
MaxTime: maximum time in seconds to perform optimization
CalcHists: at the end of the run how to output histograms: 0 to do nothing, 
           1 to compute histograms anew, 2 from reweighing
HistsIter: save histograms and potential figures every this iterations, 
           use a negative number to append iter to plot filename
"""
            
        #set initial time and other parameters
        self.CheckReady()
        StartTime = time.time()
        self.StartIter = self.Iter
        self.Mode = "INIT"
        self.Backtracking = False
        self.UseMaxChange = True
        self.PrevLineSearch = False
        self.Initdx()
        
        if self.Verbose:
            print self.Output0()
        self.OutputTarHistFile()
                           
        #initialize any constraints and reset parameters
        self.InitConstraints(Verbose = self.Verbose)
        self.SetParam(self.Param)

        #make a trajectory if necessary
        if not self.ReweightTar:
            self.UpdateModTraj()
            self.OutputModHistFile()
            self.OutputPlot()
            
        #store this mask
        NotFixed = np.logical_not(self.Fixed)

        #record the initial parameters
        self.InitParam = self.Param.copy()

        #setup for minimization
        if self.Verbose:
            print "Starting minimization\n"

        #compute the srel and obj derivatives
        self.CalcObj(CalcDeriv = True)
            
        #update log and screen
        self.OutputLog()
        if self.Verbose:
            print self.Output2()   
            
        #run convergence loop
        Converged = False
        LastHistsIter = self.StartIter
        while True:
            #check whether or not we should continue
            if not MaxTime is None:
                if time.time() - StartTime > MaxTime:
                    break
            if Converged:
                if MinIter is None:
                    break
                elif MinIter < self.Iter - self.StartIter:
                    break
            if not MaxIter is None:
                if MaxIter <= self.Iter - self.StartIter:
                    break
            
            #see whether we want to adjust dx this round
            Adjustdx = self.Autodx and self.Iter - self.StartIter > 1
            #store the derivatives
            LastDObj = self.DObj.copy()
            LastDDObj = self.DDObj.copy()
            #get ready to update the parameter set                
            self.PrepUpdate()
            #compute the new search dir
            self.CalcSearchDir(ForceSteepest = self.Iter - self.StartIter < SteepestIter)
            #get an estimate of dx
            self.Getdx()
            #make sure we don't change limits for any parameter
            DistLim, SoftDistLim = self.CalcDistLim()
            if self.dx >= DistLim:
                self.dx = 0.5 * DistLim
            
            #loop unil a good enough parameter set is found
            while True:

                #get new parameters                
                NewParam = self.Param + self.dx * self.Dir
                self.SetParam(NewParam)   
                
                #compute estimated parameter change for old param
                dP = self.Param - self.LastParam
                DeltaSrel1 = np.dot(dP, LastDObj) + 0.5 * np.dot(dP.T, np.dot(LastDDObj, dP)) 
                    
                #compute the srel and obj derivatives
                self.CalcObj(CalcDeriv = True)
                    
                #update dx if needed
                self.Backtracking = False
                if Adjustdx:
                    #compute estimated parameter change for new param
                    DeltaSrel2 = np.dot(dP, self.DObj) + 0.5 * np.dot(dP.T, np.dot(self.DDObj, dP))
                    #examine the difference; if two big, scale down the step size and backtrack
                    AutodxErr = abs(DeltaSrel1 - DeltaSrel2)
                    AutodxFracErr = AutodxErr / min(abs(DeltaSrel1), abs(DeltaSrel2))
                    if (AutodxErr < self.AutodxErrMin or AutodxFracErr < self.AutodxFracErrMin):
                        self.Setdx(ScaleFactor = self.dxScaleUp)
                    elif AutodxErr > self.AutodxErrMax and AutodxFracErr > self.AutodxFracErrMax:
                        self.Setdx(ScaleFactor = self.dxScaleDown)
                        self.Backtracking = True
                        
                if np.any(self.ReweightFrac < self.MinReweightFrac) and not self.Backtracking:
                    #moved too far away from the saved trajectory, make a new one
                    self.OutputLog()
                    if self.Verbose:
                        print self.Output2()
                        print "Reweight fraction reached %.4f; making new trajectory."  % np.min(self.ReweightFrac)
                    self.UpdateModTraj()
                    self.CalcObj(CalcDeriv = True)
                                        
                #update the log and summary file and terminal output
                self.Mode = self.DirModeStr()
                self.OutputLog()
                self.OutputSumFile()

                #output detailed histograms and potentials this iteration?
                if CalcHists > 0 and self.Iter - LastHistsIter >= abs(HistsIter):
                    LastHistsIter = self.Iter
                    self.OutputModHistFile()
                    self.OutputPlot(UseIter = HistsIter < 0)
                if self.Verbose:
                    print self.Output2()
                    
                #break out of loop if found good parameters
                if self.Backtracking:
                    self.SetParam(self.LastParam)
                else:
                    break
    
            #update which parameters have converged
            #compute absolute tolerance and fractional tolerances for each param
            AbsDiff = (self.Param - self.LastParam) / self.Scale
            FracDiff = AbsDiff / (np.abs(self.Param) + 1.e-300)
            Test1 = np.all(AbsDiff[NotFixed] < self.NRAbsTol)
            Test2 = np.all(FracDiff[NotFixed] < self.NRFracTol)
            #also check the gradients
            Test3 = np.all(np.abs((self.DObj*self.Scale)[NotFixed]) < self.NRGradTol)                
            #final convergence
            Converged = Test1 and Test2 and Test3
                
        if self.Verbose:
            Elapsed = (time.time() - StartTime) / 60.
            print "Elapsed time: %.1f min" % Elapsed
            if not Converged:
                print "Minimization did not converge in %d iters\n" % self.Iter
            else:
                print "Done with minimization\n"

        #output summaries
        self.OutputSumFile()
        self.Mode = "FINAL"
        self.OutputLog()
        if self.Verbose:
            print self.Output3()
            
        #output final histograms
        if CalcHists == 1:
            self.UpdateModTraj()
            self.OutputModHistFile()
            self.OutputPlot(UseIter = HistsIter < 0)
        elif CalcHists == 2 and self.Reweight:
            self.OutputModHistFile()
            self.OutputPlot(UseIter = HistsIter < 0)

        return Converged      



    def Run(self, *args, **kwargs):
        """Default optimization is conjugate gradient minimization."""
        return self.RunConjugateGradient(*args, **kwargs)
        
        

    def RunStages(self, *args, **kwargs):
        """Runs optimization in stages with different coefficients on the penalty terms.
StageCoefs:  list of coefficients to use in the penalties with each stage
other arguments etc. are the same as the default minimization method"""
        self.CheckReady()
        if len(self.Penalties) == 0:
            raise AttributeError("No penalties are defined for this optimizer.")
        #initialize the lagrange multipliers to zero
        for Pen in self.Penalties:
            Pen.InitializeOptimization()
        if not "StageCoefs" in kwargs:
            raise SyntaxError("Must specify input argument StageCoefs")
        StageCoefs = kwargs.pop("StageCoefs")
        for (i, Coef) in enumerate(StageCoefs):
            print "="*20 + "STAGE %d / %d" % (i+1, len(StageCoefs)) + "="*20 
            print "COEF = %12.4e\n" % Coef
            self.CGFracTol = StageCGFracTols[i]
            self.CGGradTol = StageCGGradTols[i]
            print "CGFracTol = %12.4e" %self.CGFracTol 
            print "CGGradTol = %12.4e\n" %self.CGGradTol
            for Pen in self.Penalties:
                Pen.Coef = Coef
                print "LAGMULT for %s = %12.4e\n" % (Pen.Name, Pen.LagMult)
            self.Run(*args, **kwargs)
            #update the lagrange multipliers
            for Pen in self.Penalties:
                Pen.UpdateLagMult() 
                
