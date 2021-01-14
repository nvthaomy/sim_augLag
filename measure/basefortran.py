#/usr/bin/env python

### Fortran routines for measurements in SIM suite.
### coded by MSS

import numpy as np

import sim.fortran as fortran

#Base code that potentials can define
#All measures should have a [DOMEASURE] somewhere
SourceBase = """
>>> inputmaps
>>> outputmaps
>>> defs
>>> externals
>>> init
>>> beforepreloop1
>>> preloop1beforepair
>>> preloop1pair
>>> preloop1afterpair
>>> beforemainloop
>>> mainloopbeforepair
>>> mainlooppair
>>> mainloopafterpair
>>> beforeangleloop
>>> angleloop
>>> beforetorsionloop
>>> torsionloop
>>> final
"""

SourceVars = """
>>> inputmaps
bool MeasureAll = MeasureAll
int StepNum = StepNum
int CycleNum = CycleNum
float Weight = Weight

>>> defs
bool UseMeasure(NMeasure)
bool DoMeasure1
bool DoMeasure2
int v1
int v2
float Val0
float Val1
int HistInd0
int HistInd1
"""


SourceVal_1 = """
Val0 = Val(0)
ValSum(0) = ValSum(0) + Val0*Weight
ValSumSq(0) = ValSumSq(0) + Val0*Val0*Weight
Count = Count + Weight
"""

SourceHist_1 = """
HistCount = HistCount + Weight
HistInd0 = int((Val0 - HistMin(0)) * HistiBin(0))
if (HistInd0 >= 0 .and. HistInd0 < %(NBin0)d) then
    Hist(HistInd0) = Hist(HistInd0) + Weight
endif
""" 

SourceVal_2a = """
Val0 = Val(0)
Val1 = Val(1)
ValSum(0) = ValSum(0) + Val0*Weight
ValSum(1) = ValSum(1) + Val1*Weight
ValSumSq(0,0) = ValSumSq(0,0) + Val0*Val0*Weight
ValSumSq(0,1) = ValSumSq(0,1) + Val0*Val1*Weight
ValSumSq(1,0) = ValSumSq(1,0) + Val1*Val0*Weight
ValSumSq(1,1) = ValSumSq(1,1) + Val1*Val1*Weight
Count = Count + Weight
"""

SourceVal_2b = """
Val(0) = Val0
Val(1) = Val1
ValSum(0) = ValSum(0) + Val0*Weight
ValSum(1) = ValSum(1) + Val1*Weight
ValSumSq(0) = ValSumSq(0) + Val0*Val0*Weight
ValSumSq(1) = ValSumSq(1) + Val1*Val1*Weight
Count = Count + Weight
"""

SourceHist_2 = """
HistCount = HistCount + Weight
HistInd0 = int((Val0 - HistMin(0)) * HistiBin(0))
HistInd1 = int((Val1 - HistMin(1)) * HistiBin(1))
if (HistInd0 >= 0 .and. HistInd0 < %(NBin0)d .and. HistInd1 >= 0 .and. HistInd1 < %(NBin1)d) then
    Hist(HistInd0,HistInd1) = Hist(HistInd0,HistInd1) + Weight
endif
"""

SourceVal_Nb = """
ValSum = ValSum + (Val)*Weight
ValSumSq = ValSumSq + (Val)*(Val)*Weight
Count = Count + Weight
"""

SourceVal_Na = """
ValSum = ValSum + (Val)*Weight
do v1 = 0, %(NVal)d - 1
    do v2 = 0, %(NVal)d - 1
        ValSumSq(v1,v2) = ValSumSq(v1,v2) + Val(v1)*Val(v2)*Weight
    enddo
enddo
Count = Count + Weight
"""


SourceUseMeasure = """
!measure %(MName)s
if (Active) then
    if (MeasureAll) then
        %(Use)s = .true.
        Val = 0.d0
    elseif (StepFreq > 0 .and. mod(StepNum, StepFreq)==0) then
        %(Use)s = .true.
        Val = 0.d0
    elseif (CycleFreq > 0 .and. mod(CycleNum, CycleFreq)==0) then
        %(Use)s = .true.
        Val = 0.d0
    else
        %(Use)s = .false.
    endif
else
    %(Use)s = .false.
endif"""

SourceUseMeasureGeneric = """
!measure %(MName)s
if (Active) then
    if (MeasureAll) then
        %(Use)s = .true.
    elseif (StepFreq > 0 .and. mod(StepNum, StepFreq)==0) then
        %(Use)s = .true.
    elseif (CycleFreq > 0 .and. mod(CycleNum, CycleFreq)==0) then
        %(Use)s = .true.
    else
        %(Use)s = .false.
    endif
else
    %(Use)s = .false.
endif"""


SourceMain = """
!##### compute initial quantities #####                      || HBoxL,iBoxL,DoMinImage
HBoxL = 0.5 * BoxL                                           || HBoxL
iBoxL = 1.d0 / max(1.d-300, BoxL)                            || iBoxL
DoMinImage = any(BoxL > 0.d0)                                || DoMinImage
LoopMode = 0                                                 || LoopMode
"""


def NewFortCode(Source):
    """Returns a new FortCode class."""    
    fc = fortran.FortCode(Source = Source)
    fcBase = fortran.FortCode(Source = SourceBase)
    fc.AddMissing(fcBase)
    return fc


def ReplaceTokens(Source, M, MInd):
    """Replaces all variables in source code with master arrays."""
    #replace references to measure object
    Source = fortran.ReplaceToken(Source, "THISMEASURE", "self[%d]" % MInd)
    #update do measure section
    Source = fortran.ReplaceToken(Source, "[DOMEASURE]", "[DOMEASURE%d]" % MInd)
    #update variables with module names
    for (Var, ModuleVar) in M.LibVars0.items():
        Source = fortran.ReplaceToken(Source, Var, ModuleVar)
    return Source


def AddDoMeasure(M, MInd, Md, Mfc):
    if M.IsAvg:
        if M.NVal == 1:
            s = SourceVal_1
            if hasattr(M, "Hist"):
                s += "\n" + SourceHist_1
        elif M.NVal == 2:
            if M.UseCovar:
                s = SourceVal_2a
            else:
                s = SourceVal_2b 
            if hasattr(M, "Hist"):
                s += "\n" + SourceHist_2
        else:
            if M.UseCovar:
                s = SourceVal_Na
            else:
                s = SourceVal_Nb
        s = s % Md
        s = ReplaceTokens(s, M, MInd)
        Mfc.Add(s, Block = "DOMEASURE%d" % MInd)
    return Mfc


def GetFortCode(MList, Sys):
    """Returns a FortCode object for the measurement function."""
    #fortran code object
    fc = fortran.FortCode(Source = SourceMain)
    fc.Add(fortran.loops.GetSourceVars())
    fc.Add(SourceVars)

    #make measurement items
    MGroup = []
    for (MInd,M) in enumerate(MList):
        #first check 
        fc.CheckFortCode(M.FortCode, ErrOrigin = M.Name)
        #now replace variables and make new
        Source = M.FortCode.RawSource()
        Source = ReplaceTokens(Source, M, MInd)
        Mfc = fortran.FortCode(Source, Comment = "measure %s" % M.Name, VarIndex = MInd)
        fc.AddFortCodeVars(Mfc)
        Md = {"MName":M.Name, "MInd":MInd, "NVal":getattr(M, "NVal", 0), "Use":"UseMeasure(%d)" % MInd}
        if hasattr(M, "Hist"):
            for j in range(M.NVal):
                Md["NBin%d" % j] = M.HistNBin[0]
        Mfc = AddDoMeasure(M, MInd, Md, Mfc)
        MGroup.append((M, Md, Mfc))
        
    #compute whether or not to use each measure
    for (MInd, (M, Md, Mfc)) in enumerate(MGroup):
        if M.IsAvg:
            s = SourceUseMeasure % Md
        else:
            s = SourceUseMeasureGeneric % Md
        s = ReplaceTokens(s, M, MInd)
        fc.Add(s)
    
    fc.Break()
        
    #add the init commands
    for (MInd, (M, Md, Mfc)) in enumerate(MGroup):
        s = Mfc.ProcessMain("[init]")
        s = fortran.IfBracket(s, Md["Use"])
        fc.Add(s)
        fc.Break()
                    
    #make loops
    for LoopPrefix in ["preloop1", "mainloop"]:
        
        #add the next commands
        for (MInd, (M, Md, Mfc)) in enumerate(MGroup):
            s = Mfc.ProcessMain("[before%s]" % LoopPrefix)
            s = fortran.IfBracket(s, Md["Use"])
            fc.Add(s)
       
        #make loop 
        l = fortran.loops.Loop1or2()
    
        #add commands from each poential
        UsesLoop = set()
        for (MInd, (M, Md, Mfc)) in enumerate(MGroup):
            if M.Filter is None:
                FilterN = 0
            else:
                FilterN = M.Filter.N
            s = Mfc.ProcessMain("[%sbeforepair]" % LoopPrefix)
            s = fortran.IfBracket(s, Md["Use"])
            if s: UsesLoop.add(MInd)
            if FilterN == 1: 
                l.AddLoop(SourceSingle = s, Filter = M.Filter, Sys = Sys)
            elif FilterN > 1:
                l.AddLoop(SourceSingle = s, Filter = M.Filter[0], Sys = Sys)
            s = Mfc.ProcessMain("[%spair]" % LoopPrefix)
            s = fortran.IfBracket(s, Md["Use"])
            if s: UsesLoop.add(MInd)
            if FilterN == 2:
                l.AddLoop(SourcePair = s, Filter = M.Filter, Sys = Sys)
            s = Mfc.ProcessMain("[%safterpair]" % LoopPrefix)
            s = fortran.IfBracket(s, Md["Use"])
            if s: UsesLoop.add(MInd)
            if FilterN == 1: 
                l.AddLoop(SourceAfterPair = s, Filter = M.Filter, Sys = Sys)
            elif FilterN > 1:
                l.AddLoop(SourceAfterPair = s, Filter = M.Filter[0], Sys = Sys)
        
        #add loop
        s = " .or. ".join([MGroup[i][1]["Use"] for i in sorted(UsesLoop)])
        if s: l.IfCondition = s
        fc.AddFortCode(l)
    
    #add the next commands
    for (MInd, (M, Md, Mfc)) in enumerate(MGroup):
        s = Mfc.ProcessMain("[beforeangleloop]")
        s = fortran.IfBracket(s, Md["Use"])
        fc.Add(s)
        fc.Break()
        
    #angle loop
    UsesLoop = set()
    l = fortran.loops.LoopNIntra(N=3)
    for (MInd, (M, Md, Mfc)) in enumerate(MGroup):
        if M.Filter is None: continue
        if not (M.Filter.N == 3 and M.Filter.Intra): continue
        s = Mfc.ProcessMain("[angleloop]")
        s = fortran.IfBracket(s, Md["Use"])
        l.AddLoop(s, Filter = M.Filter, Sys = Sys)
        if s: UsesLoop.add(MInd)
    
    #add angle loop
    s = " .or. ".join([MGroup[i][1]["Use"] for i in sorted(UsesLoop)])
    if s: l.IfCondition = s    
    fc.AddFortCode(l)
    
    #add the next commands
    for (MInd, (M, Md, Mfc)) in enumerate(MGroup):
        s = Mfc.ProcessMain("[beforetorsionloop]")
        s = fortran.IfBracket(s, Md["Use"])
        fc.Add(s)
        fc.Break()
    
    #torsion loop
    UsesLoop = set()
    l = fortran.loops.LoopNIntra(N=4)
    for (MInd, (M, Md, Mfc)) in enumerate(MGroup):
        if M.Filter is None: continue
        if not (M.Filter.N == 4 and M.Filter.Intra): continue
        s = Mfc.ProcessMain("[torsionloop]")
        s = fortran.IfBracket(s, Md["Use"])
        l.AddLoop(s, Filter = M.Filter, Sys = Sys)
        if s: UsesLoop.add(MInd)
    
    #add torsion loop
    s = " .or. ".join([MGroup[i][1]["Use"] for i in sorted(UsesLoop)])
    if s: l.IfCondition = s      
    fc.AddFortCode(l)   

    #add the final commands
    for (MInd, (M, Md, Mfc)) in enumerate(MGroup):
        s = Mfc.ProcessMain("[final]")
        s = fortran.IfBracket(s, Md["Use"])
        fc.Add(s)   
        fc.Break()
        
    return fc
        
