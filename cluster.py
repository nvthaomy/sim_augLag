#!/usr/bin/env python

import geom
import numpy as np
import copy


def BasicRMSD(Pos1, Pos2, CalcInd = None):
    """Calculates the rmsd between two configurations without alignment.
    * Pos1: array of dimensions [N,3] for conformation 1
    * Pos2: array of dimensions [N,3] for conformation 2
    * CalcInd: indices in [0,N) for positions in Pos to compute RMSD
"""
    if CalcInd is None:
        rmsdsq = np.sum((Pos1-Pos2)**2) / float(Pos1.shape[0])
    else:
        rmsdsq = np.sum((Pos1[CalcInd]-Pos2[CalcInd])**2) / float(len(CalcInd))
    rmsdsq = np.max([rmsdsq,0.])
    return np.sqrt(rmsdsq)
    
    
def RMSD(Pos1, Pos2, AlignPos2 = False, CompInd = None, CalcInd = None):
    """Calculates the rmsd between two configurations with alignment."
    * Pos1: array of dimensions [N,3] for conformation 1
    * Pos2: array of dimensions [N,3] for conformation 2
    * AlignPos2: True = modify Pos2 to be aligned to Pos1 (default is False)
    * CompInd: indices in [0,N) for positions in Pos to perform alignment
    * CalcInd: indices in [0,N) for positions in Pos to compute RMSD
"""
    #get indices
    if CompInd is None:
        p1, p2, n = Pos1, Pos2, len(Pos1)
    else:
        p1, p2, n = Pos1.take(CompInd, axis=0), Pos2.take(CompInd, axis=0), len(CompInd)
    #check for correct shapes
    if not p1.shape == p2.shape:
        raise ValueError("Position vectors are not the same size.")
    elif not len(p1.shape) == 2:
        raise ValueError("Position vectors are not the correct rank.")
    #get alignment
    Pos1Vec, Pos2Vec, RotMat, Resid = geom.AlignmentRMSD(p1, p2, Center = True)
    #compute rmsd
    if not CompInd is None or not CalcInd is None:
        if CalcInd is None:
            p1, p2, n = Pos1, Pos2, len(Pos1)
        else:
            p1, p2, n = Pos1.take(CalcInd, axis=0), Pos2.take(CalcInd, axis=0), len(CalcInd)
        p1 = p1 + Pos1Vec
        p2 = np.dot(p2 + Pos2Vec, RotMat)
        Resid = np.sum((p1 - p2)**2, axis=None)
    r = np.sqrt(Resid / float(n))
    #align Pos2 to Pos1
    if AlignPos2: 
        Pos2[:,:] = np.dot(Pos2 + Pos2Vec, RotMat) - Pos1Vec
    #return rmsd
    return r
    
    
def AlignPos(Pos1Ref, Pos2Ref, Pos2):
    """Aligns Pos2 using RMSD transformation matrix to match Pos2Ref to Pos1Ref."""
    #get alignment
    Pos1Vec, Pos2Vec, RotMat, Resid = geom.AlignmentRMSD(Pos1Ref, Pos2Ref, Center = True)
    return np.dot(Pos2 + Pos2Vec, RotMat) - Pos1Vec
    
    

def __MaskPos(Pos, ClustMask):
    if ClustMask is None:
        return Pos
    else:
        return Pos[ClustMask]


def ClusterMSS(Trj, Cutoff, ClustMask = None, ReportFrames = False, 
    CompInd = None, CalcInd = None, Weights = None, Verbose = True,
    MaxIter = 10, MaxCluster = None, MaxClusterWork = None, 
    IterMaxCluster = False, IterNormalize = False, Method = 0):
    """Clusters conformations in a trajectory based on RMSD distance.
Input:
    * Trj: a trajectory object
    * Cutoff: maximum RMSD distance of a configuration to a cluster
    * ClustMask: can be used to supply a mask for clustering purposes
    * ReportFrames: true to report individual stats of all frames to clusters
    * CompInd: indices in [0,N) for positions in Pos to perform alignment
    * CalcInd: indices in [0,N) for positions in Pos to compute RMSD
    * Weights: weighting factor for each conformation
    * MaxIter: maximum number of iterations to perform
    * MaxCluster: maximum number of clusters; negative values will force
      all configs to be members of a cluster (default is none)
    * MaxClusterWork: maximum number of working clusters
    * IterMaxCluster: True will dump all but MaxCluster configs each iter
    * IterNormalize: True will dump previous iter contribs to centroids
    * Method: 0 to do full rotation/alignment for RMSD; 1 to skip
Returns a tuple of:
    * Pos: list of position arrays for cluster configurations
    * ClustNum: gives the cluster number of each frame in Trj (0 for none, negative if centroid)
    * ClustWeights: gives the total weight in each cluster (fractional)
    * ClustPop: gives number of conformations in each cluster
    * PosRmsd: gives frame rmsds to corresponding clusters
    * ClustRmsd: gives cluster by cluster rmsds
    * ClustFrame: gives the frame number of trj corresponding to cluster centroids
    * ReportStr: a string reporting results
"""
    #initialize variables
    Iteration = 0      #iteration number
    WeightSum = []     #total weights of clusters
    PosSum = []        #list of cluster configuration arrays
    FinalIters = 0     #number of iterations without additions/deletions of clusters
    NTrj = len(Trj)
    if Weights is None: Weights = np.ones(NTrj, float)
    Weights = np.array(Weights, float)
    if not len(Weights) == NTrj:
        raise IndexError, "Incorrect number of array elements in Weights."
        
    #filter weights for too low values
    Weights = Weights.copy()
    Weights = Weights / Weights.max()
    Weights[Weights < 1.e-100] = 0.
    NTrjUse = int(np.sum(Weights > 0))
    
    StartInd, NewStartInd = 0, -1
    if Verbose:
        print "Using %d conformations in %d trajectory frames." % (NTrjUse, NTrj)
        
    #loop over iterations
    while FinalIters < 2:
        
        Iteration += 1
        FinalIters += 1
        if Iteration > MaxIter:
            if Verbose: print "Did not converge within maximum number of iterations"
            break
            
        if Verbose: print "Cluster iteration %d" % Iteration
        if Verbose: print "Starting with %d clusters" % len(PosSum)

        ClustNum = np.zeros(NTrj, int)  #cluster number of each configuration, starting at 1
        NAddThis = [0]*len(PosSum)      #number of configs added to each cluster this iteration
        ThisFrame = 0
        PosSumThis = copy.deepcopy(PosSum)
        WeightSumThis = copy.deepcopy(WeightSum)
        
        #check where to start
        if NewStartInd >= 0: StartInd = NewStartInd
        NewStartInd = -1
        
        #loop over configurations
        for CurInd in range(StartInd, NTrj) + range(0, StartInd):
            
            CurPos = __MaskPos(Trj.Get(CurInd), ClustMask)
            
            #check for zero weight
            CurWeight = Weights[CurInd]
            if CurWeight == 0.:
                ClustNum[CurInd] = 0
                continue
                
            ThisFrame += 1
            ind = -1  #cluster number assigned to this config; -1 means none
            
            #calculate the rmsd between this configuration and each cluster config,
            #but stop when a rmsd is found which is below the cutoff
            minRMSD = 1.e300
            for (i, PosSumi) in enumerate(PosSum):
                ThisPos = PosSumi / WeightSum[i]
                if Method == 0:
                    #rmsd between new config and average cluster config
                    r = RMSD(ThisPos, CurPos, AlignPos2 = True, 
                             CompInd = CompInd, CalcInd = CalcInd)
                else:
                    #rmsd between new config and average cluster config, without alignment
                    r = BasicRMSD(ThisPos, CurPos, CalcInd = CalcInd)
                minRMSD = min(minRMSD, r)
                if r < Cutoff:
                    #go with a cluster if rmsd is within the cutoff
                    ind = i
                    break
                    
            if ind >= 0:
                #add the configuration to the cluster
                PosSum[ind] = PosSum[ind] + CurPos * CurWeight
                WeightSum[ind] = WeightSum[ind] + CurWeight
                NAddThis[ind] = NAddThis[ind] + 1
                ClustNum[CurInd] = ind + 1
                
            elif len(PosSum) < MaxClusterWork or MaxClusterWork is None:
                #create a new cluster with this config, as long as it
                #doesn't exceed the maximum number of working clusters
                if minRMSD == 1.e300: minRMSD = 0.
                if Verbose: 
                    print "Adding cluster: config %d (%d/%d) | min RMSD %.1f | %d clusters tot" % (Trj.Index+1,
                          ThisFrame, NTrjUse, minRMSD, len(PosSum)+1)
                PosSum.append(CurPos * CurWeight)
                WeightSum.append(CurWeight)
                NAddThis.append(1)
                ClustNum[CurInd] = len(PosSum)
                FinalIters = 0
                
            else:
                #cluster is nothing
                ClustNum[CurInd] = 0
                FinalIters = 0
                if NewStartInd < 0:
                    NewStartInd = CurInd
                    if Verbose: 
                        print "Ran out of clusters. Next iteration starting from config %d" % (Trj.Index+1,)
                        
        #remove contribution to centroids from all but this round
        if IterNormalize:
            for i in range(len(PosSumThis)):
                PosSum[i] = PosSum[i] - PosSumThis[i]
                WeightSum[i] = WeightSum[i] - WeightSumThis[i]
        del PosSumThis
        del WeightSumThis
        
        #loop through clusters
        i = 0
        while i < len(PosSum):
            #remove clusters that have no additions this iteration
            if NAddThis[i] == 0:
                if Verbose: print "Removing cluster %d" % (i+1,)
                del PosSum[i]
                del WeightSum[i]
                del NAddThis[i]
                for (k, cn) in enumerate(ClustNum):
                    if cn > i + 1:
                        ClustNum[k] -= 1
                    elif cn == i + 1:
                        ClustNum[k] = -1
                FinalIters = 0
            else:
                i += 1
                
        #sort clusters and then remove any beyond MaxCluster
        ClustNum = np.array(ClustNum, int)
        PosSum, ClustNum, WeightSum, NAddThis = __SortClust(PosSum, ClustNum, Weights, WeightSum, NAddThis, Verbose)
        
        #crop off any extraneous clusters; clusterless configs
        #are assigned a cluster index of 0
        if IterMaxCluster and not MaxCluster is None and len(PosSum) > abs(MaxCluster):
            del PosSum[abs(MaxCluster):]
            WeightSum = WeightSum[:abs(MaxCluster)]
            NAddThis = NAddThis[:abs(MaxCluster)]
            ClustNum[abs(ClustNum) > abs(MaxCluster)] = 0
            
    #crop off any extraneous clusters; clusterless configs
    #are assigned a cluster index of 0
    if not IterMaxCluster and not MaxCluster is None and len(PosSum) > abs(MaxCluster):
        del PosSum[abs(MaxCluster):]
        WeightSum = WeightSum[:abs(MaxCluster)]
        NAddThis = NAddThis[:abs(MaxCluster)]
        ClustNum[abs(ClustNum) > abs(MaxCluster)] = 0
        
    #finalize things
    if Verbose: print "Calculating average structures"
    Pos = [x / y for (x, y) in zip(PosSum, WeightSum)]
    del PosSum
    del WeightSum
    
    #get cluster populations
    ClustWeights, ClustPop = __CalcClustPop(Pos, ClustNum, Weights)
    
    #if there is a maximum cluster specification that's negative, force
    #everything to the closest cluster
    if not MaxCluster == None and MaxCluster < 0:
        ClustWeights, ClustPop, Pos, ClustNum = __CalcForceClust(Trj, ClustWeights, ClustPop, Pos, ClustNum, Weights,
                                                   ClustMask, CompInd, CalcInd, Verbose)
                                                   
    #normalize weights
    ClustWeights = ClustWeights / np.sum(Weights)
    
    #calculate final rmsd values for configs and clusters
    Pos, ClustNum, PosRmsd, ClustRmsd, ClustFrame = __CalcRmsd(Trj, Pos, ClustNum, 
                                                               ClustMask, CompInd, CalcInd, Verbose)
                                                               
    #aligns second and higher cluster positions to first one
    if Method == 0:
        MaskPos0 = __MaskPos(Pos[0], ClustMask)
        NewPos = [Pos[0]] 
        for Posi in Pos[1:]:
            MaskPosi = __MaskPos(Posi, ClustMask)    
            NewPosi = AlignPos(MaskPos0, MaskPosi, Posi)
            NewPos.append(NewPosi)     
        Pos = NewPos                      
    
    #get string of results
    s = ClustResultStr(Trj, Pos, ClustNum, ClustWeights, ClustPop, 
                       PosRmsd, ClustRmsd, ClustFrame, ReportFrames)
    
    #finalize and return
    if Verbose: print "%d configurations sorted into %d clusters" % (NTrjUse, len(Pos))
    return Pos, ClustNum, ClustWeights, ClustPop, PosRmsd, ClustRmsd, ClustFrame, s


def __SortClust(Pos, ClustNum, Weights, WeightSum, NAddThis, Verbose = True):
    if Verbose: print "Reordering clusters by population"
    Sums = [(np.sum(Weights[abs(ClustNum) == i+1]), i) for i in range(len(Pos))]
    Sums.sort()
    Sums.reverse()
    #resort the arrays
    Pos = [Pos[j] for (i,j) in Sums]
    WeightSum = [WeightSum[j] for (i,j) in Sums]
    NAddThis = [NAddThis[j] for (i,j) in Sums]
    #create a dictionary which will tell us the new cluster
    #number for a given old cluster number
    Trans = {0:0}
    for i in range(len(Pos)):
        ind = Sums[i][1] + 1
        Trans[ind] = i + 1
        Trans[-ind] = -i - 1
    #update ClustNum with the rearranged cluster numbers
    ClustNum = np.array([Trans[i] for i in ClustNum], int)
    return Pos, ClustNum, WeightSum, NAddThis


def __CalcClustPop(Pos, ClustNum, Weights):
    #update the cluster population
    ClustWeights = np.array([np.sum(Weights[abs(ClustNum) == i+1]) for i in range(len(Pos))], float)
    #update the populations
    ClustPop = [float(np.sum(abs(ClustNum) == i+1)) for i in range(len(Pos))]
    ClustPop = np.array(ClustPop, float)
    return ClustWeights, ClustPop


def __CalcForceClust(Trj, ClustWeights, ClustPop, Pos, ClustNum, Weights,
                     ClustMask = None, CompInd = None, CalcInd = None, Verbose = True):
    #count the number of clusterless configurations
    c = np.sum(ClustNum == 0)
    if Verbose: print "Forcing %d extraneous configurations to existing clusters" % c
    #find the nearest cluster to each clusterless config and assign it
    for (j, CurPos) in enumerate(Trj):
        if ClustNum[j] == 0:
            ind = -1
            minr = 0.
            CurPos = __MaskPos(CurPos, ClustMask)
            for (i, Posi) in enumerate(Pos):
                r = RMSD(Posi, CurPos, AlignPos2 = False, 
                         CompInd = CompInd, CalcInd = CalcInd)
                if r < minr or ind < 0:
                    ind = i
                    minr = r
            ClustNum[j] = ind + 1
            ClustWeights[ind] = ClustWeights[ind] + Weights[j]
            ClustPop[ind] = ClustPop[ind] + 1.
    return ClustWeights, ClustPop, Pos, ClustNum


def __CalcRmsd(Trj, Pos, ClustNum,
               ClustMask = None, CompInd = None, CalcInd = None, Verbose = True):
    if Verbose: print "Calculating cluster rmsd values"
    #calculate the pairwise cluster rmsd values
    ClustRmsd = np.zeros((len(Pos), len(Pos)),float)
    for (i, Posi) in enumerate(Pos):
        for (j, Posj) in enumerate(Pos):
            if j <= i: continue
            ClustRmsd[i,j] = RMSD(Posi, Posj, AlignPos2=False, 
                                  CompInd = CompInd, CalcInd = CalcInd)
            ClustRmsd[j,i] = ClustRmsd[i,j]
    if Verbose: print "Calculating final rmsd values"
    #loop through configs and find the one with the lowest
    #rmsd in each cluster
    PosRmsd = -1. * np.ones(len(ClustNum), float)
    ClustFrame = [-1]*len(Pos)
    for (CurInd, CurPos) in enumerate(Trj):
        i = abs(ClustNum[CurInd]) - 1
        if i >= 0:
            CurPos = __MaskPos(CurPos, ClustMask)
            PosRmsd[CurInd] = RMSD(Pos[i], CurPos, AlignPos2=False, 
                                   CompInd = CompInd, CalcInd = CalcInd)
            if ClustFrame[i] < 0:
                ClustFrame[i] = CurInd
            elif PosRmsd[ClustFrame[i]] > PosRmsd[CurInd]:
                ClustFrame[i] = CurInd
    #loop through the configs again and extract the
    #coords of the minimum-rmsd configs for each clust
    if Verbose: print "Finding nearest cluster structures"
    for (i, ind) in enumerate(ClustFrame):
        ClustNum[ind] = -ClustNum[ind]
        Pos[i] = Trj.Get(ind)
    return Pos, ClustNum, PosRmsd, ClustRmsd, ClustFrame
    
    


def ClustResultStr(Trj, Pos, ClustNum, ClustWeights, ClustPop, PosRmsd, 
                   ClustRmsd, ClustFrame, ReportFrames = False):
    #make the indices
    ConfIndices = Trj.GetIndices()
    #calculate the percent
    ClustPct = [100.*w for w in ClustWeights]
    #make a results string
    s = "CLUSTER POPULATION:\nCluster number, Number of configs, Percent\n"
    s += "\n".join(["%-5d  %-7d  %.2f" % (i+1, ClustPop[i], ClustPct[i])
                    for i in range(0,len(Pos))])
    s += "\n\nCLUSTER-TO-CLUSTER RMSD:\nCluster number, Cluster number, RMSD\n"
    for i in range(0,len(Pos)):
        for j in range(i+1,len(Pos)):
            s += "%-5d  %-5d  %-8.2f\n" % (i+1, j+1, ClustRmsd[i,j])
    s += "\n\nCLUSTER CONFIG NUMBERS:\nCluster number, Config number, RMSD\n"
    for (cn, i) in enumerate(ClustFrame):
        s += "%-5d  %-5d  %-8.2f\n" % (cn+1, ConfIndices[i], PosRmsd[i])
    if ReportFrames:
        s += "\n\nCLUSTER MEMBERSHIP:\nConfig number, Cluster number, RMSD\n"
        for (i, r) in enumerate(PosRmsd):
            if r < 0.:
                s += "%-7d  %-5s  %-8s\n" % (ConfIndices[i], "--", "--")
            else:
                s += "%-7d  %-5d  %-8.2f\n" % (ConfIndices[i], ClustNum[i], PosRmsd[i])
    return s


def WriteClustResults(Trj, ClustRet, SummaryFilename, PdbFilename):
    """Writes cluster results to a text file and trajectory."""
    from sim.traj import PdbWrite
    Pos, ClustNum, ClustWeights, ClustPop, PosRmsd, ClustRmsd, ClustFrame, s = ClustRet      
    file(SummaryFilename,"w").write(s)
    pdbtrj = PdbWrite(PdbFilename)
    dummy = Trj[0]
    HeadVars = Trj.GetHeadVars()
    FrameData = Trj.FrameData.copy()
    for (ThisPos, ThisWeight) in zip(Pos, ClustWeights):
        FrameData["REMARK"] = "Cluster percentage %.2f" % (100. * ThisWeight)
        pdbtrj.WriteFrame(ThisPos, HeadVars, FrameData)
    pdbtrj.Close()
