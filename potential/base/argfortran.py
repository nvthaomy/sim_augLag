#/usr/bin/env python

### Fortran parsing / code generation for
### computing values of arguments in energy function.


import sim.fortran as fortran
import potentialtypes as ptypes

import basefortran


#Base code that potentials can define
SourceBase = """
>>> arginit
>>> argbeforepreloop
>>> argpreloopbeforepair
>>> argprelooppair
>>> argpreloopafterpair
>>> argbeforemainloop
>>> argmainloopbeforepair
>>> argmainlooppair
>>> argmainloopafterpair
>>> arglocaldensityloop1
>>> arglocaldensityloop2
>>> argbeforeangleloop
>>> argangleloop
>>> argbeforetorsionloop
>>> argtorsionloop
>>> argfinal
>>> argevalinit
>>> argevalconst
>>> argeval
"""

SourceVarsEval = basefortran.SourceVars + """
>>> defs
int ArgType
float ArgVal
int BinInd
float ThisHist
"""

SourceVarsGet = basefortran.SourceVars + """
>>> inputmaps
float Weight = Weight

>>> defs
int ArgType
float ArgVal
int BinInd
float ThisHist
"""

SourceArgGetStats = """
>>> ARGGET
if (Weight > 0.d0) then
    ArgMin(ArgType) = min(ArgMin(ArgType), ArgVal)
    ArgMax(ArgType) = max(ArgMax(ArgType), ArgVal)
    ArgCount(ArgType) = ArgCount(ArgType) + Weight
    ArgSum(ArgType) = ArgSum(ArgType) + ArgVal * Weight
    ArgSumSq(ArgType) = ArgSumSq(ArgType) + ArgVal * ArgVal * Weight 
endif
>>> ARGWEIGHT
ArgWeightSumStats = ArgWeightSumStats + Weight
"""

SourceArgGetHist = """
>>> ARGGET
if (%(HistNBin)d > 0) then
    BinInd = int((ArgVal - ArgHistMin(ArgType)) * ArgHistiBinw(ArgType))
    if (BinInd >= 0 .and. BinInd < %(HistNBin)d) ArgHist(ArgType,BinInd) = ArgHist(ArgType,BinInd) + Weight
endif
>>> ARGWEIGHT
ArgWeightSumHist = ArgWeightSumHist + Weight
"""

SourceUpdates = """
>>> UPDATEENERGY
PEnergy = PEnergy + ThisU * ThisHist
ThisPEnergy = ThisPEnergy + ThisU * ThisHist

>>> UPDATEVIRIAL
Virial = Virial + ThisW * ThisHist
"""

SourceArgEval = """
do ArgType = 0, %(ArgNType)d - 1
    [argevalconst]
    do i = 0, %(HistNBin)d - 1
        ThisHist = ArgHist(ArgType, i) 
        if (ThisHist == 0.d0) cycle
        ArgVal = ArgHistMin(ArgType) + ArgHistBinw(ArgType) * (0.5d0 + i)
        [argeval]
    enddo
enddo
"""

SourceMainGet = """
Scale = 1.d0
LoopMode = 0
"""

SourceMainEval = """
!compute initial quantities
HBoxL = 0.5 * BoxL                                                 || HBoxL
iBoxL = 1.d0 / max(1.d-300, BoxL)                                  || iBoxL
PEnergy = 0.d0
Virial = 0.d0
if (CalcDUParam) then
    DUParam = 0.d0
    DDUParam = 0.d0
endif
if (CalcDWParam) then
    DWParam = 0.d0
    DDWParam = 0.d0
endif
Terms = 0.d0
Scale = 1.d0
"""


def ReplaceTokens(Source, FF, P, PInd = None, Coef = None):
    """Replaces ArgArray variables with Fortran variables."""
    #first do base
    Source = basefortran.ReplaceTokens(Source, FF, P, PInd, Coef)
    #update variables with module vars
    for (Var, ModuleVar) in P.Arg.LibVars0.items():
        Source = fortran.ReplaceToken(Source, "Arg" + Var, ModuleVar)
    #term energy
    if not PInd is None:
        Source = fortran.ReplaceToken(Source, "ThisPEnergy", "Terms(%d)" % PInd)        
    return Source


def GetFortCodeGet(FF, SourceArgGet):
    """Returns defs and fortran code for the potential energy function."""
    #start new fortran code
    Source = SourceMainGet + fortran.loops.GetSourceVars() + SourceVarsGet
    #sub in number of potential energy terms, which remains constant
    Source = fortran.ReplaceToken(Source, "NTerm", "%d" % len(FF))
    fc = fortran.FortCode(Source = Source)
    
    #make force field code objects
    PGroup = []
    PfcTemplate = fortran.FortCode(SourceBase)
    for (PInd,P) in enumerate(FF):
        Pd = {"PName":str(P), "PInd":PInd, "HistNBin":P.Arg.HistNBin, 
              "ArgNType":P.Arg.NType, "CutSq":"CutSq(%d)" % PInd}
        #first check code
        fc.CheckFortCode(P.FortCode, ErrOrigin = P.Label)
        #now replace parameters and make new
        Pfc = P.FortCode.Copy()
        Pfc.Add(SourceUpdates)
        Pfc.Add(SourceArgGet)
        Source = Pfc.RawSource() % Pd
        Source = ReplaceTokens(Source, FF, P = P, PInd = PInd)
        Pfc = fortran.FortCode(Source, Comment = "potential %s" % P.Label)
        Pfc.AddMissing(PfcTemplate)
        Pfc.VarIndex = PInd
        fc.AddFortCodeVars(Pfc)
        PGroup.append((P, Pd, Pfc))
    
    #add the init commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[arginit]")
        fc.Add(s)
        
    #add the weight tabulating commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[ARGWEIGHT]")
        fc.Add(s)        
                    
    #make loops
    for Prefix in ["preloop", "mainloop"]:
        
        #add the next commands
        for (P,Pd,Pfc) in PGroup:
            s = Pfc.ProcessMain("[argbefore%s]" % Prefix)
            fc.Add(s)
        
        #make loop 
        l = fortran.loops.Loop1or2()
    
        #add commands from each poential
        for (P,Pd,Pfc) in PGroup:
            s = Pfc.ProcessMain("[arg%sbeforepair]" % Prefix)
            if P.Filter.N == 1: 
                l.AddLoop(SourceSingle = s, Filter = P.Filter, Sys = FF.Sys)
            elif P.Filter.N > 1:
                l.AddLoop(SourceSingle = s, Filter = P.Filter[0], Sys = FF.Sys)
            s = Pfc.ProcessMain("[arg%spair]" % Prefix)
            if P.Filter.N == 2:
                l.AddLoop(SourcePair = s, Filter = P.Filter, Sys = FF.Sys, PairCutSqStr = Pd["CutSq"])
            s = Pfc.ProcessMain("[arg%safterpair]" % Prefix)
            if P.Filter.N == 1: 
                l.AddLoop(SourceAfterPair = s, Filter = P.Filter, Sys = FF.Sys)
            elif P.Filter.N > 1:
                l.AddLoop(SourceAfterPair = s, Filter = P.Filter[0], Sys = FF.Sys)
        
        #add loop
        fc.AddFortCode(l)
        
    #make local density energy loop; can't process targets because it's multibody
    l = fortran.loops.Loop1or2(IsAdditive = False)
    for (P,Pd,Pfc) in PGroup:
        if not P.Type == ptypes.LocalDensityPotential: continue
        s = Pfc.ProcessMain("[arglocaldensityloop1]")
        l.AddLoop(SourceSingle = s, Filter = P.Filter[0], Sys = FF.Sys)
    fc.AddFortCode(l)
    
    #make local density force loop; can't process targets because it's multibody
    l = fortran.loops.Loop1or2(IsAdditive = False)
    for (P,Pd,Pfc) in PGroup:
        if not P.Type == ptypes.LocalDensityPotential: continue
        s = Pfc.ProcessMain("[arglocaldensityloop2]")
        l.AddLoop(SourcePair = s, Filter = P.Filter, Sys = FF.Sys, PairCutSqStr = Pd["CutSq"])
    fc.AddFortCode(l)
    
    #add the next commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[argbeforeangleloop]")
        fc.Add(s)     
        
    #angle loop
    l = fortran.loops.LoopNIntra(N=3)
    for (P,Pd,Pfc) in PGroup:
        if not P.Type == ptypes.AnglePotential: continue
        s = Pfc.ProcessMain("[argangleloop]")
        l.AddLoop(s, Filter = P.Filter, Sys = FF.Sys)
    
    #add angle loop
    fc.AddFortCode(l)
    
    #add the next commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[argbeforetorsionloop]")
        fc.Add(s)
    
    #torsion loop
    l = fortran.loops.LoopNIntra(N=4)
    for (P,Pd,Pfc) in PGroup:
        if not P.Type == ptypes.TorsionPotential: continue
        s = Pfc.ProcessMain("[argtorsionloop]")
        l.AddLoop(s, Filter = P.Filter, Sys = FF.Sys)
    
    #add torsion loop
    fc.AddFortCode(l)   
            
    #add the final commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[argfinal]")
        fc.Add(s)

    return fc 
    

def GetFortCodeStats(FF):
    return GetFortCodeGet(FF, SourceArgGetStats)

def GetFortCodeHist(FF):
    return GetFortCodeGet(FF, SourceArgGetHist)



def GetFortCodeEval(FF):
    """Returns defs and fortran code for calculating arguments."""
        #start new fortran code
    Source = SourceMainEval + fortran.loops.GetSourceVars() + SourceVarsEval
    #sub in number of potential energy terms, which remains constant
    Source = fortran.ReplaceToken(Source, "NTerm", "%d" % len(FF))
    fc = fortran.FortCode(Source = Source)
    
    #make force field code objects
    PGroup = []
    PfcTemplate = fortran.FortCode(SourceBase)
    for (PInd,P) in enumerate(FF):
        Pd = {"PName":str(P), "PInd":PInd, "HistNBin":P.Arg.HistNBin, 
              "ArgNType":P.Arg.NType, "CutSq":"CutSq(%d)" % PInd}
        #first check code
        fc.CheckFortCode(P.FortCode, ErrOrigin = P.Label)
        #now replace parameters and make new
        Pfc = P.FortCode.Copy()
        Pfc.Add(SourceUpdates)
        Pfc.Add(SourceArgEval)
        Source = Pfc.RawSource() % Pd
        Source = ReplaceTokens(Source, FF, P = P, PInd = PInd, Coef = "ThisHist")
        Pfc = fortran.FortCode(Source, Comment = "potential %s" % P.Label)
        Pfc.AddMissing(PfcTemplate)
        Pfc.AddMapParent("self[%d]" % PInd)
        Pfc.VarIndex = PInd
        fc.AddFortCodeVars(Pfc)
        PGroup.append((P, Pd, Pfc))
    
    fc.Break()
    
    #add the init commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[argevalinit]")
        fc.Add(s)
        fc.Break()           
        
    #add individual potentials
    for (P,Pd,Pfc) in PGroup:
        if P.Cut is None or P.Type == ptypes.LocalDensityPotential:
            SourceCut = ">>> cuttest\n"
        else:
            SourceCut = ">>> cuttest\nif (ArgVal*ArgVal > %(CutSq)s) cycle" % Pd
        s = Pfc.ProcessMain(SourceArgEval % Pd, SourceCut) 
        s = ReplaceTokens(s, FF, P, Coef = "ThisHist")       
        fc.Add(s)
        fc.Break()  
        
    return fc 
    