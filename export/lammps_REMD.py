#/usr/bin/env python

#####################################################################
####### LAMMPS EXPORT FOR REPLICA EXCHANGE MOLECULAR DYNAMICS #######
####### ADAPTED FROM SCOTT P. CARMICHAEL                      #######
#######                                                       #######
#######                                                       #######
#####################################################################

import os, time, subprocess, uuid, warnings
import numpy as np
import sim.traj
from sim.export import lammps as LammpsExport

# REMD parameters
TEMPS = None
NCores = 8; HighTemp = 500
NStepsSwap = 1000
NStepsSwapMin = 100


def pickleTraj(TrajName, LogFile = None, LogFileToken = None, Verbose = False):
    import os, cPickle as pickle
    import sim
    
    pickleName = TrajName + '.pickle'
    if os.path.isfile(pickleName):
        of = open(pickleName, 'r')
        if Verbose: print 'Loading from pickle...'
        Trj = pickle.load(of)
        of.close()
    else:
        if Verbose: print 'Pickling trajectory...'
        Trj = sim.traj.lammps.Lammps(TrajName, LogFile = LogFile, LogFileToken = LogFileToken)
        of = open(pickleName, 'w')
        pickle.dump(Trj, of)
        of.close()
    
    init = Trj[0] # needed to enable framedata parsing
    return Trj
    

def MakeLammpsReplicaMD(Sys, NStepsMin = 10000, NStepsEquil = 100000, NStepsProd = 5000000,
                 WriteFreq = 1000, Prefix = "", TrajFile = "lammpstrj",
                 LammpsCommandsBefore = "", LammpsCommandsAfter = "", *args, **kwargs):
    """Makes a LAMMPS MD trajectory, with default settings.
Returns InFile, DataFile, TableFile, DihedralFile, TrajFile"""
    
    global TEMPS, NCores, HighTemp, NStepsSwap
    import sim.units
    
    #add a prefix
    TrajFile = Prefix + TrajFile

    #start input file
    s = ""
    if len(LammpsCommandsBefore):
        s += LammpsCommandsBefore + "\n\n"

    #choose the timestep
    if Sys.Units == sim.units.DimensionlessUnits:
        TimeStep = Sys.Int.Method.TimeStep
    elif Sys.Units == sim.units.MKSUnits:
        TimeStep = Sys.Int.Method.TimeStep
    elif Sys.Units == sim.units.AtomicUnits:
        TimeStep = Sys.Int.Method.TimeStep * 1000. / 20.5  #convert to fs
    else:
        raise LammpsError("Don't recognize type of units in system.")

    #make a dictionary for filling in template
    RanSeed = 100000 * np.random.rand()
    RanSeed2 = 100000 * np.random.rand()
    RanSeed3 = 100000 * np.random.rand()

    #use exponentially distributed temperatures
    if TEMPS is None: TEMPS = np.logspace(np.log10(Sys.TempSet),np.log10(HighTemp), NCores)
    else: NCores = len(TEMPS)

    #routine checks on swap freq
    if NStepsSwap < WriteFreq:
        print "Note: Swapping faster than writing co-ordinates to file"
    if NStepsSwap < NStepsSwapMin:
        print "WARNING: Swapping too frequently. Increase NStepsSwap"
     
    d = {"TIMESTEP" : TimeStep,
         "TEMP" : Sys.TempSet,
         "LANGEVINDAMP" : 1. / Sys.Int.Methods.VVIntegrate.LangevinGamma,
         "WRITEFREQ" : WriteFreq,
         "NSTEPSMIN" : NStepsMin,
         "NSTEPSEQUIL" : NStepsEquil,
         "NSTEPSPROD" : NStepsProd,
         "NSTEPSTEMPER" : NStepsEquil + NStepsProd,
         "NSTEPSSWAP" : NStepsSwap,
         "TRAJFILE" : TrajFile,
         "RANSEED" : RanSeed,
         "RANSEED2" : RanSeed2,
         "RANSEED3" : RanSeed3,
         "TEMPS" : ' '.join( [ '%11.4f'% i for i in TEMPS ] ),
         "REPINDS": ' '.join( str(i) for i in range( len(TEMPS) ) )
         }

    #add to input string
    s += """
#temperature schedule
variable        t world %(TEMPS)s
variable        rep world %(REPINDS)s

#minimization
minimize        1.e-4 0.0 %(NSTEPSMIN)d %(NSTEPSMIN)d

#timestep
timestep        %(TIMESTEP)-11.4e

#time integration
fix             timeintegration all nve

#thermostat
fix             thermostat all langevin $t $t %(LANGEVINDAMP)11.4e %(RANSEED)d
""" %d

    #recentering if not periodic
    if np.any(Sys.BoxL <= 0):
        s += "fix             recenterfix all recenter 0.0 0.0 0.0 units box\n"

        #run equilibration and production (REMD)
    s += """
#restart thermo
reset_timestep  0
thermo          %(WRITEFREQ)d
thermo_style    custom etotal ke pe temp ebond eangle edihed evdwl
thermo_modify   line multi norm no format float %%%%14.7e

#setup trajectory output
dump            myDump all custom %(WRITEFREQ)d %(TRAJFILE)s.${rep}.gz id type x y z element
dump_modify     myDump element ATOMNAMES
dump_modify     myDump sort id
dump_modify     myDump format  "%%%%d %%%%d %%%%14.7f %%%%14.7f %%%%14.7f %%%%s"

#run production
temper          %(NSTEPSTEMPER)d %(NSTEPSSWAP)d $t thermostat %(RANSEED2)d %(RANSEED3)d
""" % d

    if len(LammpsCommandsAfter):
        s += "\n\n" + LammpsCommandsAfter + "\n"

    #now make the lammps files
    LammpsFiles = LammpsExport.MakeLammps(Sys, *args, Prefix = Prefix, LammpsCommands = s, **kwargs)
    #return filenames
    return LammpsFiles, TrajFile


def RunLammpsReplica(InFile, Prefix = "", LogFile = "lammps.log", Verbose = False,
              CheckOutput = True):
    """Runs REMD using lammps.  Returns LogFile, ScreenFile, ReturnCode."""

    global TEMPS, NCores
    slaveID = Prefix + uuid.uuid4().hex[:8]
    LogFile = Prefix + LogFile
    ScreenFile = Prefix + 'screen'

    d = {'SLAVEID': slaveID, 'JOBNAME': slaveID.split('/')[-1],
         'PREFIX': Prefix,
         'SLAVEDIR': os.getcwd(),
         'SLAVEFILE': os.path.join(os.path.dirname(InFile), slaveID + '.sh'),
         'HEADNODE': 'zin.cnsi.ucsb.edu',
         'NCORES': NCores,
         'LMPEXEC': LammpsExport.LammpsExec,
         'INFILE': InFile, 'LOGFILE': LogFile, 'SCREENFILE': ScreenFile}

    s_Slave = """
#!/bin/bash
#
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -N %(JOBNAME)s
#$ -pe mpich %(NCORES)d
mpirun -np %(NCORES)d %(LMPEXEC)s -partition %(NCORES)dx1 -in %(INFILE)s -log %(LOGFILE)s -screen %(SCREENFILE)s
"""
    s_Monitor = """
# mark the start of the run
touch .started
# if running from head node, submit job straightaway
# else export all relevant SGE variables from current compute node to head node and submit from there
if [ $HOSTNAME = %(HEADNODE)s ]; then OUT=$(qsub %(SLAVEFILE)s);
else OUT=$(ssh %(HEADNODE)s "export SGE_ROOT=/opt/gridengine ; export SGE_QMASTER_PORT=536 ; cd %(SLAVEDIR)s ; $SGE_ROOT/bin/linux-x64/qsub %(SLAVEFILE)s")
fi

# set up a monitor until job is done
JOBID=$(echo $OUT | awk '{print $3}')
while qstat -u $USER | grep $JOBID &> /dev/null; do
  sleep 2;
done
if [ ! -d .trash ]; then mkdir .trash; fi
mv %(SLAVEFILE)s .trash/
mv %(SLAVEID)s.o$JOBID .trash/
mv %(SLAVEID)s.po$JOBID .trash/
mv .started .trash/
"""
    file(d['SLAVEFILE'], 'w').write(s_Slave % d) 
    if Verbose: print "Starting LAMMPS REMD job..."
    t1 = time.time()
    p = subprocess.Popen(s_Monitor % d, shell = True, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    p.communicate()
    returncode = p.returncode #TODO: get this returncode from Lammps
    if Verbose:
        t2 = time.time()
        print "LAMMPS finished in %d seconds." % (t2-t1)
    if CheckOutput: LammpsExport.CheckLammpsOutput(LogFile)
    return LogFile, ScreenFile, returncode


def MakeLammpsReplicaTraj(Sys, DelTempFiles = False, Prefix = "", Verbose = False, *args, **kwargs):
    """Returns a LAMMPS trajectory object.  Takes same arguments as MakeLammpsReplicaMD.
Returns Traj, TrajFile."""

    global TEMPS, NCores, HighTemp, NStepsSwap

    RepInds = []
    LammpsFiles, TrajFile = MakeLammpsReplicaMD(Sys, Prefix = Prefix, *args, **kwargs)
    InFile, DataFile, TableFile, DihedralFile = LammpsFiles
    LogFile, ScreenFile, returncode = RunLammpsReplica(InFile, Prefix = Prefix, Verbose = Verbose)

    NStepsEquil = kwargs.get('NStepsEquil', 0)
    NStepsProd = kwargs.get('NStepsProd', 0)
    WriteFreq = kwargs.get('WriteFreq', 1)

    if Verbose: print 'Reordering trajectories at %g' % Sys.Temp
    TempInd = [list(TEMPS).index(T) for T in list(TEMPS) if '%3.2f' % T  == '%3.2f' % Sys.TempSet][0]
    RepIndsMaster = [np.where(x[1:] == TempInd)[0][0] for x in np.loadtxt(LogFile, skiprows = 3)]
    RepInds = RepIndsMaster[int(NStepsEquil/NStepsSwap) : ]
    RepInds = RepInds[:-1]               
    this_Traj = {} ; TrajList = [] ; MultiTrajFnList = []
    
    if Verbose: t1 = time.time()
    if WriteFreq <= NStepsSwap: # assume mod(NStepsSwap, StepFreq) = 0
        # extract only the prod run replica indices
        for ii, i in enumerate(RepInds):
            if not i in this_Traj.keys():
                thisTrajFn = '%s.%d.gz' % (TrajFile, i)
                thisLogFn = '%s.%d' % (LogFile, i)
                this_Traj[i] = pickleTraj(thisTrajFn, LogFile = thisLogFn, LogFileToken = '#run production')
            # get last NStepsProd of this_Traj to write to temp. ordered traj
            this_Traj_slice = this_Traj[i][int(NStepsEquil/WriteFreq) : ]
            start = ii * NStepsSwap / WriteFreq
            stop = (ii + 1) * NStepsSwap / WriteFreq
            TrajList.append(this_Traj_slice[start:stop])
    
    else:
        NSkip = WriteFreq / NStepsSwap # assume mod(StepFreq, NStepsSwap) = 0
        for ii, i in enumerate(RepInds[0::NSkip]):
            if not i in this_Traj.keys():
                thisTrajFn = '%s.%d.gz' % (TrajFile, i)
                thisLogFn = '%s.%d' % (LogFile, i)
                this_Traj[i] = pickleTraj(thisTrajFn, LogFile = thisLogFn, LogFileToken = '#run production')
            this_Traj_slice = this_Traj[i][int(NStepsEquil/WriteFreq) : ]
            TrajList.append(this_Traj_slice[ii:ii+1])

    MultiTraj = sim.traj.Multi(TrajList, Sys)
    print "<<<traj length,",len(MultiTraj)
    MultiTrajFnList = ['%s.%d.gz' % (TrajFile, i) for i in range(len(TEMPS))]
    if Verbose:
        t2 = time.time()
        print 'Reordering finished in %d seconds' % (t2-t1) 
    
    if DelTempFiles:
        for fn in InFile, DataFile, TableFile, DihedralFile, LogFile, ScreenFile:
            if os.path.isfile(fn): os.remove(fn)

        for i in range(len(TEMPS)):
            this_logfile = '%s.%d' % (LogFile, i)
            this_screenfile = '%s.%d' % (ScreenFile, i)
            for fn in this_logfile, this_screenfile:
                if os.path.isfile(fn): os.remove(fn)

    return MultiTraj, MultiTrajFnList
