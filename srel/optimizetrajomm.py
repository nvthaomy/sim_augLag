#/usr/bin/env python


### Relative entropy optimizer for saved trajectories that uses OpenMM to
### make model trajectories.
### coded by KS


import optimizetraj
import base

OpenMMDelTempFiles = False
OpenMMStepsMin = 0

#1 for errors, 2 for warnings, 0 for nothing
OpenMMErrorDiffEne = 2


def MakeModTrajOpenMM(ModSys, FilePrefix, StepsEquil, StepsProd, ElecSys, StepsStride = 1, Verbose = False):
    """Makes a trajectory for the model system using OpenMM and returns ModTraj.
Note: ModTraj must expose a FileName attribute and the Delete() method."""
    #import accessory modules
    from sim.export import omm
    from sim.export import lammps
    import sim.traj as traj
    import os
    #run openmm
    if ElecSys != None: #running trajectory on ElecSys MN 2020.01.05
        print('Making CG trajectory based on a pseudo system: {}'.format(ElecSys.Name))
        r = omm.MakeOpenMMTraj(ElecSys, DelTempFiles = OpenMMDelTempFiles,
                              Prefix = FilePrefix + ".", Verbose = Verbose,
                              NStepsMin = OpenMMStepsMin, NStepsEquil = StepsEquil,
                              NStepsProd = StepsProd, WriteFreq = StepsStride)
        Traj, ElecTrajFile, ElecTrajFileDCD = r
        #reprocess trajectory and get potential energy with forcefield of ModSys 
        print('\nReprocess trajectory using forcefield of {}'.format(ModSys.Name))
        Traj = omm.MakeRerunOpenMMTraj(ModSys, ElecTrajFile, ElecTrajFileDCD, DelTempFiles = OpenMMDelTempFiles,
                              Prefix = FilePrefix + "Rerun.", Verbose = Verbose,
                              NStepsMin = OpenMMStepsMin, NStepsEquil = StepsEquil,
                              NStepsProd = StepsProd, WriteFreq = StepsStride)
        TrajFn = ElecTrajFile
    else:
        r = omm.MakeOpenMMTraj(ModSys, DelTempFiles = OpenMMDelTempFiles, 
                              Prefix = FilePrefix + ".", Verbose = Verbose,
                              NStepsMin = OpenMMStepsMin, NStepsEquil = StepsEquil,
                              NStepsProd = StepsProd, WriteFreq = StepsStride)
        #=== Post-processing to get lammpstrj that sim can read===
        Traj, TrajFn, TrajFnDCD = r

    #get energies from trajectory
    ret = base.GetEneTraj(Traj, ModSys, ErrorNoEne = 1, ErrorDiffEne = OpenMMErrorDiffEne)
    TrajEne, SysEne, TrajLogWeight, Nullx = ret
    
    #convert to raw format for faster reading/parsing
    if Verbose:
        print "Converting OpenMM trajectory to raw format for faster access."
    ModTrajFn = FilePrefix + ".traj.dat"
    traj.Convert(Traj, traj.RawWrite, ModTrajFn, FrameDataFields = ["BoxL"], 
                 NewFrameData = {"PEnergy" : SysEne}, Verbose = Verbose)
    del Traj  
#    os.remove(TrajFn) 
    return traj.Raw(ModTrajFn) 


def UseOpenMM(Opt):
    """Modifies OptimizeTrajClass instance Opt to use OpenMM."""
    if not isinstance(Opt, optimizetraj.OptimizeTrajClass):
        raise TypeError("Opt must be an OptimizeTrajClass class or subclass.")
    Opt.MakeModTrajFunc = MakeModTrajOpenMM
    return Opt
    
    
def OptimizeTrajOpenMMClass(*args, **kwargs):
    """Wrapper to produce a OpenMM-enabled OptimizeTrajClass."""
    Opt = optimizetraj.OptimizeTrajClass(*args, **kwargs)
    return UseOpenMM(Opt)
    
    
    

