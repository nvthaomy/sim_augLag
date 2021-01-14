#/usr/bin/env python


### Relative entropy optimizer for saved trajectories that uses LAMMPS to
### make model trajectories.
### coded by MSS


import optimizetraj
import base

LammpsDelTempFiles = True
LammpsStepsMin = 100

#1 for errors, 2 for warnings, 0 for nothing
LammpsErrorDiffEne = 2


def MakeModTrajLammps(ModSys, FilePrefix, StepsEquil, StepsProd, ElecSys, StepsStride = 1, Verbose = False):
    """Makes a trajectory for the model system using LAMMPS and returns ModTraj.
Note: ModTraj must expose a FileName attribute and the Delete() method."""
    #import accessory modules
    from sim.export import lammps
    import sim.traj as traj
    import os
    if ElecSys != None: #running trajectory on ElecSys MN 2020.01.05
        print('Making CG trajectory based on a pseudo system: {}'.format(ElecSys.Name))
        r = lammps.MakeLammpsTraj(ElecSys, DelTempFiles = LammpsDelTempFiles,
                              Prefix = FilePrefix + ".", Verbose = Verbose,
                              NStepsMin = LammpsStepsMin, NStepsEquil = StepsEquil,
                              NStepsProd = StepsProd, WriteFreq = StepsStride)
        Traj, ElecTrajFile = r
        #reprocess trajectory and get potential energy with forcefield of ModSys
        Traj = lammps.MakeRerunLammpsTraj(ModSys, ElecTrajFile, DelTempFiles = LammpsDelTempFiles,
                              Prefix = FilePrefix + "Rerun.", Verbose = Verbose,
                              NStepsMin = LammpsStepsMin, NStepsEquil = StepsEquil,
                              NStepsProd = StepsProd, WriteFreq = StepsStride)
        TrajFn = ElecTrajFile
    else:
        #run lammps
        #turn on Pressure
        r = lammps.MakeLammpsTraj(ModSys, DelTempFiles = LammpsDelTempFiles, 
                              Prefix = FilePrefix + ".", Verbose = Verbose,
                              NStepsMin = LammpsStepsMin, NStepsEquil = StepsEquil,
                              NStepsProd = StepsProd, WriteFreq = StepsStride)
        Traj, TrajFn = r
    #get energies from trajectory
    ret = base.GetEneTraj(Traj, ModSys, ErrorNoEne = 1, ErrorDiffEne = LammpsErrorDiffEne)
    TrajEne, SysEne, TrajLogWeight, Nullx = ret
    #convert to raw format for faster reading/parsing
    if Verbose:
        print "Converting LAMMPS trajectory to raw format for faster access."
    ModTrajFn = FilePrefix + ".traj.dat"
    traj.Convert(Traj, traj.RawWrite, ModTrajFn, FrameDataFields = ["BoxL"], 
                 NewFrameData = {"PEnergy" : SysEne}, Verbose = Verbose)
    del Traj  
    os.remove(TrajFn) 
    return traj.Raw(ModTrajFn) 


def UseLammps(Opt):
    """Modifies OptimizeTrajClass instance Opt to use LAMMPS."""
    if not isinstance(Opt, optimizetraj.OptimizeTrajClass):
        raise TypeError("Opt must be an OptimizeTrajClass class or subclass.")
    Opt.MakeModTrajFunc = MakeModTrajLammps
    return Opt
    
    
def OptimizeTrajLammpsClass(*args, **kwargs):
    """Wrapper to produce a LAMMPS-enabled OptimizeTrajClass."""
    Opt = optimizetraj.OptimizeTrajClass(*args, **kwargs)
    return UseLammps(Opt)
    
    
    

