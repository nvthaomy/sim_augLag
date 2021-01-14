#/usr/bin/env python


### Velocity verlet MD integration routine objects in SIM suite.
### coded by MSS

import numpy as np
import sim.fortran as fortran

import base.methodbase as methodbase
import velverletfortran as vvfortran


  

class VVIntegrate(methodbase.IntegrateMethodClass):
    
    ModuleVars = ["TimeStep", "RattleMaxIterPerAtom", "RattleTol1", "RattleTol2", "RattleTol3",
                  "Thermostat", "AndersenStep", "AndersenStepFreq", "AndersenCollisionFreq",
                  "Glogs", "Vlogs", "Xlogs", "QMass", "LangevinGamma", 
                  "Barostat", "BarostatStepFreq", "BarostatDeltaV", "BarostatUseAxis", 
                  "BarostatDoIsotropic", "BarostatMaxVol", "BarostatMinVol", "BarostatNAtt", "BarostatNAcc",
                  "RemoveCOMStep", "RemoveCOMStepFreq", "TEnergySum", "TEnergySqSum"]
    ModuleConsts = ["NH_N"]

    def __init__(self, Sys, TimeStep = None):
        """Initializes a velocity verlet MD integrator."""
        object.__setattr__(self, "LibVars", {})
        methodbase.IntegrateMethodClass.__init__(self, Sys, StepsPerCycle = 10)
        #fortran subroutine for integrating and steepest descent
        PreCallCode = "Sys = self.Sys"
        self.Sub0 = fortran.Subroutine("vvupdatekenergy", SubStr = "self.Sub0",
                                       PreCallCode = PreCallCode)
        self.Sub0.PreLoad = self.PreLoad
        self.Sub0.PostLoad = self.PostLoad   
        self.Sys.Lib.append(self.Sub0)
        self.Sub1 = fortran.Subroutine("vvintegrate", SubStr = "self.Sub1",
                                       PreCallCode = PreCallCode)
        self.Sys.Lib.append(self.Sub1)
        #set the timestep
        self.TimeStep = TimeStep
        if self.TimeStep == None:
            self.TimeStep = Sys.Units.TimeScale * 0.001
        #rattle tolerances and maximum number of rattle iters
        self.RattleMaxIterPerAtom = 100
        self.RattleTol1 = 1.e-4
        self.RattleTol2 = 1.e-4
        self.RattleTol3 = 1.e-4
        #add thermostat codes
        for (n,v) in vvfortran.Thermostats.iteritems():
            setattr(self, n, v)
        #default thermostat
        self.Thermostat = self.ThermostatNone
        #add barostat codes
        for (n,v) in vvfortran.Barostats.iteritems():
            setattr(self, n, v)
        #default barostat
        self.Barostat = self.BarostatNone  
        #andersen parameters
        self.AndersenStep = 0
        self.AndersenStepFreq = 100
        self.AndersenCollisionFreq = 0.1
        #setup default collision frequency
        self.SetupAndersenCollisionFreq()
        #initialize Nose Hoover
        self.SetupNoseHoover(NThermostat = 2)
        #setup default langevin gamma
        self.LangevinGamma = 0.0002 / self.TimeStep
        #removing center of mass 
        self.RemoveCOMStep = 0
        self.RemoveCOMStepFreq = 10000
        self.RemoveCOMOnInit = True
        #total energy check
        self.TEnergySum = 0.
        self.TEnergySqSum = 0.
        #default barostat frequency
        self.BarostatStepFreq = 100
        #barostat delta V default
        self.BarostatDeltaV = 10. * Sys.Units.Angstrom**self.Sys.Dim
        #axes to allow in barostat
        self.BarostatUseAxis = np.ones(self.Sys.Dim, dtype=bool)
        #isotropic barostat moves?
        self.BarostatDoIsotropic = True
        #barostat min and max volumes
        self.BarostatMaxVol = 1.e300
        self.BarostatMinVol = 0.
        #attempts and acceptances for barostat (computed each cycle)
        self.BarostatNAtt = 0.
        self.BarostatNAcc = 0.
        self.BarostatAccRatio = 0.
        
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
        

    def SetupAndersenCollisionFreq(self, Coef = 0.15, N = None, TimeStep = None):
        """Estimates/sets the value of the Andersen per particle collision rate according to
the equation nu = Coef / (TimeStep * N^(2/3)).  Estimates are made
for each parameter from current Sys conditions if not given.  Technically this
formula is valid only in 3D."""
        if TimeStep is None:
            TimeStep = self.TimeStep
        if N is None:
            if self.Sys.NDOF == 0:
                N = 200
            else:
                N = self.Sys.NDOF / self.Sys.Dim
        self.AndersenCollisionFreq = Coef / (TimeStep * N**(2./3.))
        
        
    def SetupNoseHoover(self, NThermostat = 1, Masses = None):
        """Initializes Nose-Hoover routines."""
        self.NH_N = NThermostat
        self.Glogs = np.ones(NThermostat, float)
        self.Vlogs = np.ones(NThermostat, float)
        self.Xlogs = np.ones(NThermostat, float)
        self.QMass = np.ones(NThermostat, float)
        if Masses is None:
            w = 0.01 / self.TimeStep
            self.QMass.fill(self.Sys.Units.EScale / w**2)
        elif type(Masses) is float:
            self.QMass.fill(Masses)
        elif len(Masses) == len(self.QMass):
            self.QMass[:] = Masses
        else:
            raise ValueError("Cannot process Masses; need %d values." % len(self.QMass))
        #check masses
        if np.any(self.QMass <= 0.):
            raise ValueError("One of more Nose Hoover masses is <= 0.")

    def Reset(self):
        self.AndersenStep = 0
        self.RemoveCOMStep = 0

    def Cleanup(self):
        methodbase.IntegrateMethodClass.Cleanup(self)
        self.Sub0 = None
        self.Sub1 = None
  
    def PreLoad(self):
        """Run before fortran compilation."""
        #initialize the main variables           
        self.Sys.Lib.VarPreLoad(self, self.ModuleVars, self.ModuleConsts, Prefix = "VV_")        
        #main source for energy routine
        fc0 = vvfortran.GetKEFortCode(self.Sys, self.LibVars0)
        self.Sub0.AddFortCode(fc0)
        fc1 = vvfortran.GetVVIFortCode(self.Sys, self.LibVars0)
        self.Sub1.AddFortCode(fc1)
        
    def PostLoad(self):
        """Run after fortran compilation."""
        self.Sys.Lib.VarPostLoad(self)  

    def UpdateKEnergy(self):
        """Updates kinetic energy for system."""
        exec(self.Sub0.CallObj)

    def RunAndersenThermostat(self):
        """Redraws all velocities to system, with zero COM."""
        mask = []
        for (ind, Mol) in enumerate(self.Sys):
            mask.extend([self.Sys.MolActive[ind] > 0] * len(Mol))
        mask = np.array(mask)
        p = self.Sys.Vel * self.Sys.Mass[:,np.newaxis]
        Sigma = np.sqrt(self.Sys.TempSet * self.Sys.Units.kB) * self.Sys.sqrtMass[mask]
        p[mask] = np.random.normal(size = (np.sum(mask), self.Sys.Dim)) * Sigma[:,np.newaxis]      
        p = p - p.mean(axis=0)
        self.Sys.Vel = p * self.Sys.iMass[:,np.newaxis]
        #update kinetic energy
        exec(self.Sub0.CallObj)
        
    def RunRemoveCOM(self):
        """Adjusts velocities to remove any center of mass motion."""
        #calculate the momentum 
        p = self.Sys.Mass[:,np.newaxis] * self.Sys.Vel
        #subtract off mean
        p = p - p.mean(axis=0)
        #set velocities
        self.Sys.Vel = p * self.Sys.iMass[:,np.newaxis]
        #update kinetic energy
        exec(self.Sub0.CallObj)

    def Init(self, NSteps):
        """Prepares to take MD steps."""
        self.Sys.Flags.ConservesMomentum = vvfortran.ThermostatConservesMomentum[self.Thermostat]
        self.Sys.Flags.CalcForce = True
        self.Sys.ForceField.Eval()
        if self.RemoveCOMOnInit:
            self.RunRemoveCOM()
        else:
            self.UpdateKEnergy()
        self.TEnergySum = 0.
        self.TEnergySqSum = 0.
        ##set the minimum volume to be the case when the minimum box length
        #is twice the biggest cutoff
        BoxL = [x for x in self.Sys.BoxL if x > 0.]
        if len(BoxL):
            minBoxL = min(BoxL)
        else:
            minBoxL = 1.
            if not self.Barostat == self.BarostatNone:
                raise ValueError("Barostat is set but no box length is supplied.")          
        maxCut = 0.
        for P in self.Sys.ForceField:
            if not P.Cut is None:
                maxCut = max(P.Cut, maxCut)
        CurVol = np.prod(self.Sys.BoxL)
        self.BarostatMinVol = 1.0000001 * CurVol * (2. * maxCut / minBoxL)**self.Sys.World.Dim
        #zero acceptance statistics
        self.BarostatNAtt = 0.
        self.BarostatNAcc = 0.
       
    def Cycle(self):
        """Completes one Velocity Verlet integration cycle."""
        self.ErrorFlag = 0
        exec(self.Sub1.CallObj)
        if self.ErrorFlag > 0:
            self.ErrorMsg = vvfortran.ErrorMsgs.get(self.ErrorFlag, "Integration error")
            return 0
        return self.StepsPerCycle

    def Finalize(self, NSteps):
        """Finalizes after MD steps."""
        self.TEnergyAvg = self.TEnergySum / float(NSteps)
        self.TEnergyStd = self.TEnergySqSum / float(NSteps) - self.TEnergyAvg**2
        self.TEnergyStd = np.sqrt(max(self.TEnergyStd, 0.))
        self.BarostatAccRatio = self.BarostatNAcc / (self.BarostatNAtt + 1.e-10)



def FindTimeStep(Sys, NSteps, EneFracError = 1.e-5, Verbose = True,
                 GuessTimeStep = None, ScaleUp = 1.2, ScaleDown = 0.5, 
                 MaxIter = 20, ErrorTol = 0.1):
    """Attempts to find an appropriate time step."""
    Int = Sys.Int
    Method = Int.Method
    if not isinstance(Method, VVIntegrate):
        raise TypeError("Can only run this function when VVIntegrate is the integrator method.")
    if not GuessTimeStep is None: Method.TimeStep = GuessTimeStep
    #first run an equilibration with Anderson
    pt = None
    if Verbose: print "Running equilibration with Nose-Hoover thermostat."
    Method.Thermostat = Method.ThermostatNoseHoover
    Method.Reset()
    if Verbose: pt = "Running NVT with time step %11.4e" % Method.TimeStep
    Int.Run(NSteps, ProgressText = pt)
    #save positions and velocities
    Pos, Vel = Sys.Pos.copy(), Sys.Vel.copy()
    #now switch to NVE and begin measuring time steps
    Method.Thermostat = Method.ThermostatNone
    dts, Errs = [], []
    while len(dts) < MaxIter:
        if len(dts) > 0:
            ErrMax = max(Errs) 
            ErrMin = min(Errs)
            dtMax = max(dts)
            dtMin = min(dts)
            if ErrMin > EneFracError:
                Method.TimeStep = dtMin * ScaleDown
            elif ErrMax < EneFracError:
                Method.TimeStep = dtMax * ScaleUp
            else:
                m, b = np.polyfit(np.log(Errs), np.log(dts), 1)
                Method.TimeStep = np.exp(m * np.log(EneFracError) + b)
        if Verbose: pt = "Running NVE with time step %11.4e" % Method.TimeStep
        Method.Reset()
        Sys.Pos, Sys.Vel = Pos, Vel
        Int.Run(NSteps, ProgressText = pt)
        if Method.TEnergyStd == 0:
            Err = 1.e-14
        else:
            Err = np.abs(Method.TEnergyStd / Method.TEnergyAvg)
        if Verbose: print "Found %11.4e frac energy error using timestep %11.4e" % (Err, Method.TimeStep)
        dts.append(Method.TimeStep)
        Errs.append(Err)
        if np.abs(Err - EneFracError) <= ErrorTol * EneFracError:
            if Verbose: print "Within tolerance ... done searching for timesteps."
            break
    if Verbose:
        l = sorted(zip(dts, Errs))
        print "Final list of timesteps and energy fractional errors:"
        for (dt, Err) in l:
            print "%11.4e %11.4e" % (dt, Err)            
        


class VVQuench(methodbase.IntegrateMethodClass):
    
    ModuleVars = ["TimeStep", "RattleMaxIterPerAtom", "RattleTol1", "RattleTol2", "RattleTol3"]
    ModuleConsts = []

    def __init__(self, Sys, TimeStep = None):
        """Initializes a velocity verlet quench integrator."""
        object.__setattr__(self, "LibVars", {})
        methodbase.IntegrateMethodClass.__init__(self, Sys, StepsPerCycle = 10)
        #fortran subroutine for integrating and steepest descent
        PreCallCode = "Sys = self.Sys\n"
        self.Sub1 = fortran.Subroutine("vvquench", SubStr = "self.Sub1", 
                                       PreCallCode = PreCallCode)
        self.Sub1.PreLoad = self.PreLoad
        self.Sub1.PostLoad = self.PostLoad
        self.Sys.Lib.append(self.Sub1)
        #set the timestep
        self.TimeStep = TimeStep
        if self.TimeStep == None:
            self.TimeStep = Sys.Units.TScale * 0.001
        #maximum number of rattle iters
        self.RattleMaxIterPerAtom = 100
        #rattle tolerances
        self.RattleTol1 = 1.e-4
        self.RattleTol2 = 1.e-4
        self.RattleTol3 = 1.e-4
        
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

    def Cleanup(self):
        methodbase.IntegrateMethodClass.Cleanup(self)
        self.Sub1 = None
  
    def PreLoad(self):
        """Run before fortran compilation."""
        #initialize the main variables          
        self.Sys.Lib.VarPreLoad(self, self.ModuleVars, self.ModuleConsts, Prefix = "VVQ_")
        #main source for argument routine
        fc1 = vvfortran.GetVVQFortCode(self.Sys, self.LibVars0)
        self.Sub1.AddFortCode(fc1)

        
    def PostLoad(self):
        """Run after fortran compilation."""
        self.Sys.Lib.VarPostLoad(self) 

    def Init(self, NSteps):
        """Prepares to take quench steps."""
        self.Sys.Flags.CalcForce = True
        self.Sys.ForceField.Eval()

    def Cycle(self):
        """Make one steepest descent (VV) integration cycle."""
        self.ErrorFlag = 0
        exec(self.Sub1.CallObj)
        if self.ErrorFlag > 0:
            self.ErrorMsg = vvfortran.ErrorMsgs.get(self.ErrorFlag, "Integration error")
            return 0
        return self.StepsPerCycle



