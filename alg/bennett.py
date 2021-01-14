#/usr/bin/env python

### Bennett free energy calculations.
### coded by MSS

import os, atexit
import numpy as np
import sim.traj as traj
import sim.utility as utility

DEBUG = False


class BennettException(Exception):
    pass


def fd(x, w):
    return 1. / (np.exp(w) + np.exp(x))
    
def logsumexp(x):
    shift = np.max(x)
    return np.log( np.sum( np.exp(x - shift) ) ) + shift
    
    
def BennettErr(BE11, BE21, BE22, BE12, FE):
    """Computes the error in the Bennett calculation, sqrt(var(ln(q2/q1))).
Based on Shirts et al, PRL 91, 140601 (2003)"""
    #compute all "work" measurements
    W = np.concatenate((BE12 - BE11, BE22 - BE21))
    #add free energy (ln (Q2/Q1))
    W = W - FE
    n = len(W)
    c = np.max(np.abs(W))
    terms = 1. / (2.*np.exp(-c) + np.exp(W-c) + np.exp(-W-c))
    err = (np.exp(c) / np.mean(terms) - 4.) / n
    err = np.sqrt(err)
    return err    
    

def BennettFE(BE11, BE21, BE22, BE12, Beta1 = 1, Beta2 = 1, LogWeight1 = None, LogWeight2 = None,
              Verbose = False, Tol = 1.e-10, MaxIter = 1000, Reweight = False):
    """Returns -log(Q2/Q1) = beta2 F2 - beta1 F1 using the Bennett method.  
BEij is the list of beta*E values for system i configs evaluated in system j.
If Reweight = True, then also returns two sets of log weights for reweighting
system 2 to system 1 conditions."""
    #check array sizes
    N1 = len(BE11)
    if len(BE12) != len(BE11):
        raise ValueError("BE11 and BE12 arrays are not the same length.")
    N2 = len(BE21)
    if len(BE22) != len(BE21):
        raise ValueError("BE21 and BE22 arrays are not the same length.")        
    if not (np.all(np.isfinite(BE11)) and np.all(np.isfinite(BE12)) 
            and np.all(np.isfinite(BE21)) and np.all(np.isfinite(BE22))):
        if DEBUG: 
            print "BE11 BE12 BE22 BE21"
            for i in xrange(N1):
                print i, BE11[i], BE12[i], BE22[i], BE21[i]
        raise ValueError("Found non-finite value in BE11, BE12, BE22, or BE21")
    #nondimensionalized if not already
    BE11 = Beta1 * BE11
    BE22 = Beta2 * BE22
    BE12 = Beta2 * BE12
    BE21 = Beta1 * BE21
    logN1 = np.log(N1)
    logN2 = np.log(N2)
    #weights
    if LogWeight1 is None: LogWeight1 = np.zeros(N1)
    LogWeight1 = LogWeight1 - logsumexp(LogWeight1)
    if LogWeight2 is None: LogWeight2 = np.zeros(N2)
    LogWeight2 = LogWeight2 - logsumexp(LogWeight2)
    #fe perturbation for initial guesses
    FE1 = -logsumexp(BE11 - BE12 + LogWeight1)  
    FE2 = logsumexp(BE22 - BE21 + LogWeight2)
    if Verbose:
        print "Bennett initial guesses: %f, %f" % (FE1, FE2)
    #setup first points: x is beta2 F2 - beta1 F1, y should be zero
    xo = FE1
    NVals = BE21 - BE22 + xo 
    DVals = BE12 - BE11 - xo 
    Nshift = -NVals.max()
    Dshift = -DVals.max()
    NVals = np.log( fd(NVals+Nshift, Nshift) ) + LogWeight2
    DVals = np.log( fd(DVals+Dshift, Dshift) ) + LogWeight1
    yo = logsumexp(NVals) - logsumexp(DVals) + Nshift - Dshift
    x = FE2
    y = 1.
    #loop until tolerance is met
    Iter = 0
    while np.abs(y) > Tol and Iter < MaxIter:
        Iter += 1
        #evaluate y for current x
        NVals = BE21 - BE22 + x 
        DVals = BE12 - BE11 - x 
        Nshift = -NVals.max()
        Dshift = -DVals.max()
        NVals = np.log( fd(NVals+Nshift, Nshift) ) + LogWeight2
        DVals = np.log( fd(DVals+Dshift, Dshift) ) + LogWeight1
        y = logsumexp(NVals) - logsumexp(DVals) + Nshift - Dshift
        #predict new x
        xn = (y*xo-yo*x)/(y-yo)
        xo = x
        x = xn
        yo = y
        #print messages
        if Verbose:
            print "Bennett iteration %d: current error is %.3e" % (Iter, np.abs(y))
        if np.isnan(y):
            raise BennettException("Bennett iterations led to NAN after %d iterations." % Iter)
    if Iter >= MaxIter:
        raise BennettException("Bennett iterations did not converge to %.3e in %d iterations." % (Tol, MaxIter))
    #now compute the estimated error
    FE = xo
    FEerr = BennettErr(BE11, BE21, BE22, BE12, FE)
    if Verbose:
        print "Bennett final free energy: %f +- %f" % (FE, FEerr)
    if not Reweight:
        return FE, FEerr
    AllEnergiesIn1 = np.concatenate((BE11, BE21))
    AllEnergiesIn2 = np.concatenate((BE12, BE22))
    #-beta(U2 - U1) + beta(A2-A1)
    ExpTerm = (AllEnergiesIn1 - AllEnergiesIn2) + FE + logN2 - logN1
    #find max term
    MaxTerm = -np.max(ExpTerm)
    Num = fd(ExpTerm + MaxTerm, MaxTerm)
    LogProb = np.log(Num)
    LogProb[:N1] = LogProb[:N1] + LogWeight1 + logN1
    LogProb[N1:] = LogProb[N1:] + LogWeight2 + logN2
    LogWeight = LogProb - logsumexp(LogProb)
    if any(np.isnan(LogWeight)):
        raise BennettException("Bennett iterations produced NAN in LogWeights.")
    return FE, FEerr, LogWeight[:N1], LogWeight[N1:]
    
    
def TestBennett(N1 = 1000000, N2 = 1000000):
    """Tests the Bennett routines using Gaussians."""
    avg1, sig1 = 0., 1.
    avg2, sig2 = 1., 2.5
    FE1 = -np.log(sig1)
    FE2 = -np.log(sig2)
    xvals1 = np.random.normal(avg1, sig1, N1)
    xvals2 = np.random.normal(avg2, sig2, N2)
    BE11 = (xvals1 - avg1)**2 / (2. * sig1**2)
    BE12 = (xvals1 - avg2)**2 / (2. * sig2**2)
    BE22 = (xvals2 - avg2)**2 / (2. * sig2**2)
    BE21 = (xvals2 - avg1)**2 / (2. * sig1**2)
    FEalg, FEerr, lw1, lw2 = BennettFE(BE11, BE21, BE22, BE12, Verbose = True, Reweight = True)
    print "TEST1" + "="*70, "\n"
    print "Actual free energy difference    : ", FE2 - FE1
    print "Calculated free energy difference: ", FEalg
    #make a new system 1 using reweighted 2
    xvals12 = np.concatenate((xvals1, xvals2))
    lw12 = np.concatenate((lw1, lw2))
    hist, edges = np.histogram(xvals12, bins=100, range=(avg1-3*sig1, avg1+3*sig1), 
                               weights=np.exp(lw12), density=True)
    hist0, edges = np.histogram(xvals1, bins=100, range=(avg1-3*sig1, avg1+3*sig1),
                                density=True)
    print "Histogram:"
    for (i,h) in enumerate(hist):
        print 0.5 * (edges[i] + edges[i+1]), h, hist0[i]
    

def CalcFEInt(Sys, StepsEquil, StepsProd, StepsStride = 10, 
              Param1 = None, Param2 = None, Temp1 = None, Temp2 = None, 
              NStatesInit = 2, MaxDeltaFE = 1.e300, MaxFEerr = 0.01, 
              BennettTol = 1.e-10, BennettMaxIter = 10000,
              Verbose = True, TempFileDir = None, CalcSrel = False):
    """Computes the free energy difference for taking the system from (Param1, Temp1) to
(Param2, Temp2).  Here, by free energy, meant is the log partition ratio log(Q2/Q1).
The reaction coordinate for the free energy path is denoted by Lambda and 
(Param, Temp) are linearly interpolated along it.  Iteratively bisects the reaction 
coordinate, computing pairwise free energy differences, until each interval's free
energy difference is less than MaxDeltaFE or error less than MaxFEErr."""
    ParamInit = Sys.ForceField.Param.Val.copy()
    
    #check for inputs
    NStatesInit = max(NStatesInit, 2)
    if Param1 is None and Param2 is None:
        Param1 = np.array(ParamInit)
        Param2 = np.array(ParamInit)
    elif Param1 is None or Param2 is None:
        raise ValueError("Both Param1 and Param2 must be specified together.")
    if len(Param1) != len(Param2) or len(Param1) != len(ParamInit):
        raise ValueError("Param1, Param2, Sys.ForceField.Param.Val are not all same size.")
    if Temp1 is None and Temp2 is None:
        Temp1 = Sys.TempSet
        Temp2 = Sys.TempSet
    elif Temp1 is None or Temp2 is None:
        raise ValueError("Both Temp1 and Temp2 must be specified together.")
        
    #initialize the list of lambdas, trajectory files, trajectories, free energies, and FE err
    States = []
    for i in range(NStatesInit):
        State = (float(i) / float(NStatesInit-1), None, None, None, None)
        States.append(State)
        
    #make a function to set the state of a system based on lambda
    def SetState(Lambda):
        Sys.TempSet = (1. - Lambda) * Temp1 + Lambda * Temp2
        Sys.ForceField.SetParam((1. - Lambda) * Param1 + Lambda * Param2)
        
    #make a function to delete the trajectory files upon exit
    def AtExitFunc():
        for (Lambda, TrjFn, Trj, FE, FEerr) in States:
            if not Trj is None and not TrjFn is None:
                Trj.Close()
                if os.path.isfile(TrjFn):
                    os.remove(TrjFn)
                    print "Deleted temporary trajectory file %s" % TrjFn
    atexit.register(AtExitFunc)
        
    #now iterate finding free energies
    while True:
        
        if Verbose:
            print "This iteration is starting with %d states." % len(States)
        
        #first make any missing trajectories
        for (i, (Lambda, TrjFn, Trj, FE, FEerr)) in enumerate(States):
            if TrjFn is None:
                #set system properties
                SetState(Lambda)
                #make a new trajectory
                TrjFn = utility.NewTempFile(Suffix = ".dat", Dir = TempFileDir)
                Trj = traj.RawWrite(TrjFn)
                States[i] = (Lambda, TrjFn, Trj, FE, FEerr)
                if Verbose:
                    print "Making new trajectory for lambda = %f as %s" % (Lambda, TrjFn)
                #equilibrate
                Sys.Int.Run(StepsEquil, ProgressText = "Equilibration")
                Trj.AddAction(Sys.Int, StepFreq = StepsStride)
                Sys.Int.Run(StepsProd, ProgressText = "Production")
                Trj.DelAction(Sys.Int)
                Trj.Close()
                Trj = traj.Raw(TrjFn)
                States[i] = (Lambda, TrjFn, Trj, FE, FEerr)
            
        #first prepare a system for computing energies
        #set all algorithm flags to none
        Sys.Flags.CalcsOff()
        #turn off measures
        Sys.Measures.Reset()
        Sys.Measures.AllOff()
        #now compute any missing free energy differences
        for i in range(len(States) - 1):
            j = i + 1
            (iLambda, iTrjFn, iTrj, iFE, iFEerr) = States[i]
            (jLambda, jTrjFn, jTrj, jFE, jFEerr) = States[j]
            iTemp = (1. - iLambda) * Temp1 + iLambda * Temp2
            jTemp = (1. - jLambda) * Temp1 + jLambda * Temp2
            iBeta = 1. / (Sys.Units.kB * iTemp)
            jBeta = 1. / (Sys.Units.kB * jTemp)
            if not iFE is None: continue
            #initialize arrays
            n = len(iTrj)
            BEii = np.zeros(n, float)
            BEij = np.zeros(n, float)
            BEji = np.zeros(n, float)
            BEjj = np.zeros(n, float)
            #loop through first trajectory and compute energies of second state
            if Verbose:
                print "Parsing trajectories for lambda1 = %f and lambda2 = %f" % (iLambda, jLambda) 
            SetState(jLambda)
            prog = utility.ProgressBar("Reading trajectory 1/2", n)
            for (k, Pos) in enumerate(iTrj):
                #get the energy for system 1 evaluated in system 1
                BEii[k] = iTrj.FrameData["PEnergy"] * iBeta
                #update the model positions 
                Sys.Pos = Pos
                #get the energy for system 1 evaluated in system 2
                Sys.ForceField.Eval()
                BEij[k] = Sys.PEnergy * jBeta
                #update the progress bar
                prog.Update(k)
                if DEBUG:
                    if np.isnan(BEii[k]) or np.isnan(BEij[k]):
                        print "Frame", k
                        print "PEnergyii", iTrj.FrameData["PEnergy"]
                        print "PEnergyij", Sys.PEnergy
                        print "Pos"
                        for i,p in enumerate(Pos):
                            print i, p
                        exit()
            prog.Clear()      
            SetState(iLambda)
            prog = utility.ProgressBar("Reading trajectory 2/2", n)
            for (k, Pos) in enumerate(jTrj):
                #get the energy for system 2 evaluated in system 2
                BEjj[k] = jTrj.FrameData["PEnergy"] * jBeta
                #update the model positions 
                Sys.Pos = Pos
                #get the energy for system 2 evaluated in system 1
                Sys.ForceField.Eval()
                BEji[k] = Sys.PEnergy * iBeta
                #update the progress bar
                prog.Update(k)
                if DEBUG:
                    if np.isnan(BEjj[k]) or np.isnan(BEji[k]):
                        print "Frame", k
                        print "PEnergyjj", jTrj.FrameData["PEnergy"]
                        print "PEnergyji", Sys.PEnergy
                        print "Pos"
                        for i,p in enumerate(Pos):
                            print i, p
                        exit()
            prog.Clear() 
            #run the bennett algorithm
            if Verbose:
                print "Running Bennett algorithm"
            iFE, iFEerr = BennettFE(BEii, BEji, BEjj, BEij, Verbose = Verbose, 
                                    Tol = BennettTol, MaxIter = BennettMaxIter)
            if Verbose:
                print "This free energy difference is %f" % iFE
            #update states data
            States[i] = (iLambda, iTrjFn, iTrj, iFE, iFEerr)
            
        #decide whether or not we need to bisect intervals
        NewStates = []
        n = len(States)
        MaxFEerr0 = MaxFEerr / np.sqrt(n-1)
        for (i, (Lambda, TrjFn, Trj, FE, FEerr)) in enumerate(States):
            if i < n - 1 and (np.abs(FE) > MaxDeltaFE or FEerr > MaxFEerr0):
                #we need to bisect this interval
                #first signal to recompute fe of this state
                NewStates.append((Lambda, TrjFn, Trj, None, None)) 
                #now add the intermediate state
                InterLambda = 0.5 * Lambda + 0.5 * States[i+1][0]
                NewStates.append((InterLambda, None, None, None, None))
                if Verbose:
                    print "Added new state with lambda = %f" % InterLambda
            else:
                #just add the old state
                NewStates.append((Lambda, TrjFn, Trj, FE, FEerr))
                
        #check to see if we can stop
        if len(NewStates) == len(States):
            if Verbose:
                print "Done bisecting states."
            break
        else:
            States = NewStates
            
    #return a list of free energy and lambda values
    if Verbose:
        print "Finished free energy calculation using %d states." % len(States)
        
    Lambdas = np.array([State[0] for State in States], float)
    DeltaFEs = np.array([State[-2] for State in States[:-1]], float) 
    DeltaFEerrs = np.array([State[-1] for State in States[:-1]], float)       
    NetFE = np.sum(DeltaFEs)
    NetFEerr = np.sqrt(np.sum(DeltaFEerrs**2))

    #see if we want to compute srel values
    if CalcSrel:
        #first prepare a system for computing energies
        #set all algorithm flags to none
        Sys.Flags.CalcsOff()
        #turn off measures
        Sys.Measures.Reset()
        Sys.Measures.AllOff()
        #compute Srel values with state 1 as target
        if Verbose:
            print "Calculating Srel values for state 1"
        (iLambda, iTrjFn, iTrj, iFE, iFEerr) = States[0]
        iTemp = (1. - iLambda) * Temp1 + iLambda * Temp2
        iBeta = 1. / (Sys.Units.kB * iTemp)
        Srel1 = []
        for j in range(len(States)):
            (jLambda, jTrjFn, jTrj, jFE, jFEerr) = States[j]
            jTemp = (1. - jLambda) * Temp1 + jLambda * Temp2
            jBeta = 1. / (Sys.Units.kB * jTemp)
            #initialize arrays
            n = len(iTrj)
            BEii = 0.
            BEij = 0.
            #loop through first trajectory and compute energies of second state
            if Verbose:
                print "Parsing trajectory for lambda = %f" % (jLambda) 
            SetState(jLambda)
            prog = utility.ProgressBar("Reading trajectory %d/%d" % (j+1, len(States)), n)
            for (k, Pos) in enumerate(iTrj):
                #get the energy for system 1 evaluated in system 1
                BEii += iTrj.FrameData["PEnergy"] * iBeta
                #update the model positions 
                Sys.Pos = Pos
                #get the energy for system 1 evaluated in system 2
                Sys.ForceField.Eval()
                BEij += Sys.PEnergy * jBeta
                #update the progress bar
                prog.Update(k)
            prog.Clear() 
            #normalize energies
            BEii = BEii / n
            BEij = BEij / n
            #compute free energy difference
            FE = np.sum(DeltaFEs[0:j])
            #update Srel array
            Srel1.append(BEij - BEii - FE)
        #compute srel values with state 2 as target
        if Verbose:
            print "Calculating Srel values for state 2"
        (iLambda, iTrjFn, iTrj, iFE, iFEerr) = States[-1]
        iTemp = (1. - iLambda) * Temp1 + iLambda * Temp2
        iBeta = 1. / (Sys.Units.kB * iTemp)
        Srel2 = []
        for j in range(len(States)):
            (jLambda, jTrjFn, jTrj, jFE, jFEerr) = States[j]
            jTemp = (1. - jLambda) * Temp1 + jLambda * Temp2
            jBeta = 1. / (Sys.Units.kB * jTemp)
            #initialize arrays
            n = len(iTrj)
            BEii = 0.
            BEij = 0.
            #loop through first trajectory and compute energies of second state
            if Verbose:
                print "Parsing trajectory for lambda = %f" % (jLambda) 
            SetState(jLambda)
            prog = utility.ProgressBar("Reading trajectory %d/%d" % (j+1, len(States)), n)
            for (k, Pos) in enumerate(iTrj):
                #get the energy for system 1 evaluated in system 1
                BEii += iTrj.FrameData["PEnergy"] * iBeta
                #update the model positions 
                Sys.Pos = Pos
                #get the energy for system 1 evaluated in system 2
                Sys.ForceField.Eval()
                BEij += Sys.PEnergy * jBeta
                #update the progress bar
                prog.Update(k)
            prog.Clear() 
            #normalize energies
            BEii = BEii / n
            BEij = BEij / n
            #compute free energy difference
            FE = -np.sum(DeltaFEs[j:])
            #update Srel array
            Srel2.append(BEij - BEii - FE)
        #convert to arrays
        Srel1 = np.array(Srel1, float)
        Srel2 = np.array(Srel2, float)
                    
    #delete temporary files    
    for (Lambda, TrjFn, Trj, FE, FEerr) in States:
        Trj.Close()
        if os.path.isfile(TrjFn):
            os.remove(TrjFn)
            print "Deleted temporary trajectory file %s" % TrjFn
            
    #return values
    if CalcSrel:
        return NetFE, NetFEerr, Lambdas, DeltaFEs, DeltaFEerrs, Srel1, Srel2
    else:
        return NetFE, NetFEerr, Lambdas, DeltaFEs, DeltaFEerrs
    


    

