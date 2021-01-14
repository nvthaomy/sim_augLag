#/usr/bin/env python


### Relative entropy optimizer for multiple trajectories.
### coded by MSS

import numpy as np
import optimize, optimizetraj



class OptimizeMultiTrajClass(optimizetraj.OptimizeTrajClass):
    
    def __init__(self, OptimizeTrajList, Weights = None, FilePrefix = None, 
                 Verbose = True, **kwargs):
        """Initializes a multiple Srel-minimization class.
OptimizeTrajList: list of OptimizeTraj classes. 
Weights: weights for contributions of each OptimizeTraj object to Srel sum
FilePrefix: file prefix for output
"""        
        optimize.OptimizeClass.__init__(self, [x.ModSys for x in OptimizeTrajList], [x.ElecSys for x in OptimizeTrajList], **kwargs)
        self.OptimizeTrajList = OptimizeTrajList
        #weights
        if Weights is None:
            self.Weights = np.ones(len(self.OptimizeTrajList), dtype=float)
        elif not len(Weights) == len(self.OptimizeTrajList):
            raise ValueError("Length of Weights must be the same as that of OptimizeTrajList")
        else:
            self.Weights = np.array(Weights, dtype=float)
        #set file names
        self.SumFileText = False
        self.FilePrefix = FilePrefix
        #set verbose options
        self.Verbose = Verbose
        #add defaults to class
        for (attr, val) in optimizetraj.Defaults.items():
            setattr(self, attr, val)
        #parsed yet?
        self.TarIsParsed = False
        #backtracking currently?
        self.Backtracking = False
        #first iteration
        self.StartIter = 0  
        #array of reweighting fractions
        self.ReweightFrac = np.ones(len(self.OptimizeTrajList), dtype = float)
        #check saving routines
        FilePrefixes = []
        for Opt in self.OptimizeTrajList:
            if not Opt.FilePrefix is None:
                if Opt.FilePrefix in FilePrefixes:
                    raise ValueError("More than one Optimizer has same FilePrefix.")
                FilePrefixes.append(Opt.FilePrefix)
                
    def ParseTarData(self):
        """Reads target trajectory, parses argument data, and sets argument mins & maxes. 
"""
        for Opt in self.OptimizeTrajList:
            if Opt.TarIsParsed:
                print "WARNING: OptimizeMultiTraj will re-parse argument ranges for each component optimizer."
                
        GlobalMins = None
        GlobalMaxs = None
        #run through component optimizers and find global arg mins and maxes
        for (i, Opt) in enumerate(self.OptimizeTrajList):
            fn = self.FilePrefix + "_m%d_argdata.pickle" % i
            Opt.FindTarArgRanges(ArgDataFileName = fn)
            Mins, Maxs = Opt.ModSys.ForceField.Arg.GetArgRanges()
            if GlobalMins is None:
                GlobalMins, GlobalMaxs = Mins, Maxs
            else:
                for (GlobalMin, OptMin) in zip(GlobalMins, Mins):
                    GlobalMin[:] = np.minimum(GlobalMin, OptMin)
                for (GlobalMax, OptMax) in zip(GlobalMaxs, Maxs):
                    GlobalMax[:] = np.maximum(GlobalMax, OptMax) 
                    
        #now run through and set each optimizer to the same ranges, followed by parsing
        from sim.potential.base.argarray import AddToData
        WeightsNorm = self.Weights / self.Weights.sum()
        AvgData = None
        for (i, Opt) in enumerate(self.OptimizeTrajList):
            fn = self.FilePrefix + "_m%d_argdata.pickle" % i
            Opt.ModSys.ForceField.Arg.SetHistRanges(GlobalMins, GlobalMaxs)        
            Opt.ParseTarData(AutoHistRanges = False, ArgDataFileName = fn)
            #get current data and add it to the total using the weights
            ThisData = Opt.ModSys.ForceField.Arg.GetData()
            AvgData = AddToData(AvgData, ThisData, WeightsNorm[i])
                           
        #set the individual histograms to all match the average
        for Opt in self.OptimizeTrajList:
            Opt.ModSys.ForceField.Arg.SetData(AvgData)  
            
        #no make a spline histogram for use in constraints
        #THIST MUST BE DONE WITH COMMON HISTOGRAM RANGES PER ABOVE
        self.SplineKnotHists = self.GetSplineKnotHists()    
        self.FLSplineKnotHists = self.GetFLSplineKnotHists()

        #ANY POTENTIAL ESTIMATES MUST BE RUN IMMEDIATELY AFTER                   
        self.TarIsParsed = True  
              
        
        
    def GetSplineKnotHists(self):
        ret = None
        for Opt in self.OptimizeTrajList:
            if ret is None:
                ret = Opt.GetSplineKnotHists()
            else:
                for (i, Hist) in enumerate(Opt.GetSplineKnotHists()):
                    if not Hist is None:
                        ret[i] += Hist
        return ret
    
    def GetFLSplineKnotHists(self):
        ret = None
        for Opt in self.OptimizeTrajList:
            if ret is None:
                ret = Opt.GetFLSplineKnotHists()
            else:
                for (i, Hist) in enumerate(Opt.GetFLSplineKnotHists()):
                    if not Hist is None:
                        ret[i] += Hist
        return ret    


    def UpdateModTraj(self):
        """Makes trajectories for each system."""
        for (i, Opt) in enumerate(self.OptimizeTrajList):
            if self.Verbose:
                print "Model trajectory #%d/%d" % (i+1, len(self.OptimizeTrajList))
            Opt.UpdateModTraj()      

    def DelModTraj(self):
        """Deletes any model trajectories."""
        for Opt in self.OptimizeTrajList:
            Opt.DelModTraj()
                
    def DelTarTraj(self):
        """Deletes any target trajectories."""
        for Opt in self.OptimizeTrajList:
            Opt.DelTarTraj()        
            
    def Cleanup(self):
        """Cleans up this class for deletion."""
        for Opt in self.OptimizeTrajList:
            Opt.Cleanup()
        
    
    def GetTarHists(self, HistLabel = "target"):
        """Gets a list of target histogram objects."""
        #collect and add together hists
        Hists = []
        for (i, Opt) in enumerate(self.OptimizeTrajList):
            l = "%s%d" % (HistLabel, i)
            for (j, Hist) in enumerate(Opt.GetTarHists(HistLabel = l)):
                if i==0:
                    Hists.append(Hist)
                else:
                    Hists[j] = Hists[j] + Hist
        return Hists
        
    def GetModHists(self, HistLabel = "model"):
        """Gets a list of target histogram objects."""
        #collect and add together hists
        Hists = []
        for (i, Opt) in enumerate(self.OptimizeTrajList):
            l = "%s%d" % (HistLabel, i)
            for (j, Hist) in enumerate(Opt.GetModHists(HistLabel = l)):
                if i==0:
                    Hists.append(Hist)
                else:
                    Hists[j] = Hists[j] + Hist
        return Hists            
            
    def GetBothHists(self, TarHistLabel = "target", ModHistLabel = "model"):
        """Returns target plus model hists."""
        Hists = []
        for (i, Opt) in enumerate(self.OptimizeTrajList):
            ml = "%s%d" % (ModHistLabel, i)
            tl = "%s%d" % (TarHistLabel, i)
            ThisHists = Opt.GetBothHists(TarHistLabel = tl, ModHistLabel = ml)
            for (j, Hist) in enumerate(ThisHists):
                if i==0:
                    Hists.append(Hist)
                else:
                    Hists[j] = Hists[j] + Hist
        return Hists            

    
    def TestObj(self, *args, **kwargs):
        print "="*20 + "TEST FOR MULTI OPTIMIZER" + "="*20
        optimizetraj.OptimizeTrajClass.TestObj(self, *args, **kwargs)
        n = len(self.OptimizeTrajList)
        for (i, Opt) in enumerate(self.OptimizeTrajList): 
            print "="*20 + "TEST FOR SUB OPTIMIZER %d/%d" % (i+1,n) + "="*20
            Opt.TestObj(*args, **kwargs)


    def OutputLog(self, Flush = True):
        """Update the log file after an iteration."""
        OtherVals = []
        for (i, Opt) in enumerate(self.OptimizeTrajList):
            OtherVals.append(("Srel%d" % i, Opt.Srel))
            OtherVals.append(("Bias%d" % i, Opt.Bias))
            OtherVals.append(("Obj%d" % i, Opt.Obj))
        optimizetraj.OptimizeTrajClass.OutputLog(self, Flush = Flush, OtherVals = OtherVals)
            

    def OutputPlot(self, UseIter = False, FilePrefix = None):
        """Updates the plot of histograms and distributions."""
        if FilePrefix is None:
            FilePrefix = self.FilePrefix
        if FilePrefix is None or not optimizetraj.EnablePlots:
            return
        import plot_forcefield
        PlotFmt = optimizetraj.PlotFmt
        #output the potentials and histograms
        if not self.ModHistFn is None:
            for (ind, Opt) in enumerate(self.OptimizeTrajList):
                if UseIter:
                    PlotFn = "%s_potentials%06d_m%d.%s" % (FilePrefix, self.Iter, ind, PlotFmt)
                else:
                    PlotFn = "%s_potentials_m%d.%s" % (FilePrefix, ind, PlotFmt)
                plot_forcefield.PlotForceField(self.ModSys, PlotFn, self.ModHistFn, ind)
            #now output the total histograms
            if UseIter:
                PlotFn = "%s_potentials%06d_tot.%s" % (FilePrefix, self.Iter, PlotFmt)
            else:
                PlotFn = "%s_potentials_tot.%s" % (FilePrefix, PlotFmt)
            plot_forcefield.PlotForceField(self.ModSys, PlotFn, self.ModHistFn, -1)                
        #output the time progression
        if not self.LogFn is None and self.Iter > self.StartIter:
            PlotFn = "%s_iterations.%s" % (FilePrefix, optimizetraj.PlotFmt)
            plot_forcefield.PlotParameterHistory(self.LogFn, PlotFn)
            

    def RevertSys(self):          
        pass 
            
            
    def CalcObj(self, CalcDeriv = True):
        """Computes the objective function and derivatives by calling member objects."""
        self.Obj = 0.
        self.DObj.fill(0.)
        self.DDObj.fill(0.)
        self.Srel = 0.
        self.Bias = 0.
        for (Opt, Weight) in zip(self.OptimizeTrajList, self.Weights):
            Opt.CalcObj(CalcDeriv = CalcDeriv)
            self.Srel += Opt.Srel * Weight
            self.Obj += Opt.Obj * Weight
            self.DObj += Opt.DObj * Weight
            self.DDObj += Opt.DDObj * Weight
            self.Bias += Opt.Bias * Weight  
            self.FluctBias += Opt.FluctBias * Weight
        #compute reweighting fraction
        self.ReweightFrac = np.array([Opt.ReweightFrac for Opt in self.OptimizeTrajList], dtype=float)
        #calculate objective and compensate for constraints
        self.UpdateObj() 
        
    
    def CheckReady(self):
        if not self.TarIsParsed:
            print "ParseTarData() not previously called.  Now calling..."
            self.ParseTarData()        
        for Opt in self.OptimizeTrajList:
            Opt.CheckReady()
        
    def RunStages(self, *args, **kwargs):
        """Runs optimization in stages with different coefficients on the penalty terms.
StageCoefs:  list of coefficients to use in the penalties with each stage
other arguments etc. are the same as the default minimization method"""
        self.CheckReady()
        HasPenalties = False
        for Opt in self.OptimizeTrajList:
            if len(Opt.Penalties) > 0:
                HasPenalties = True
        if not HasPenalties:
            raise AttributeError("No penalties found in any of the optimizer objects.")
        #initialize the lagrange multipliers to zero
        for Opt in self.OptimizeTrajList:
            for Pen in Opt.Penalties:
                Pen.InitializeOptimization()
        if not "StageCoefs" in kwargs:
            raise SyntaxError("Must specify input argument StageCoefs")
        StageCoefs = kwargs.pop("StageCoefs")
        for (i, Coef) in enumerate(StageCoefs):
            print "="*20 + "STAGE %d / %d" % (i+1, len(StageCoefs)) + "="*20 
            print "COEF = %12.4e\n" % Coef
            nopt = len(self.OptimizeTrajList)
            for (ind, Opt) in enumerate(self.OptimizeTrajList):
                print "Optimizer %d of %s" % (ind+1, nopt)
                for Pen in Opt.Penalties:
                    Pen.Coef = Coef
                    print "  LAGMULT for %s = %12.4e\n" % (Pen.Name, Pen.LagMult)
            self.Run(*args, **kwargs)
            #update the lagrange multipliers
            for Opt in self.OptimizeTrajList:
                for Pen in Opt.Penalties:
                    Pen.UpdateLagMult()         
        
