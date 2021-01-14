#/usr/bin/env python

### Fortran parsing / code generation for
### potential energy function in SIM suite.


import sim.fortran as fortran
import potentialtypes as ptypes


#Base code that potentials can define
SourceBase = """
>>> init
>>> beforepreloop
>>> preloopbeforepair
>>> prelooppair
>>> preloopafterpair
>>> beforemainloop
>>> mainloopbeforepair
>>> mainlooppair
>>> mainloopafterpair
>>> localdensityloop1
>>> localdensityloop2
>>> beforeangleloop
>>> angleloop
>>> beforetorsionloop
>>> torsionloop
>>> final
>>> saveenergystate
>>> revertenergystate
"""

SourceVars = """
>>> externals

>>> inputmaps
int Mode = Mode
bool CalcForce = Sys.Flags.CalcForce
bool CalcVirial = Sys.Flags.CalcVirial
bool CalcDUParam = Sys.Flags.CalcDUParam
bool CalcDWParam = Sys.Flags.CalcDWParam
bool CalcFluct = Sys.Flags.CalcFluct

>>> defs
float Scale
float Forcei(Dim)
float Forcek(Dim)
float Forcej(Dim)
float Forcel(Dim)
float ThisU
float ThisW
float ThisdU
float ThisA
float ThisForce(Dim)
float ThisForcesum(Dim)
"""


SourceUpdates = """
>>> UPDATEENERGY
PEnergy = PEnergy + ThisU * Scale
ThisPEnergy = ThisPEnergy + ThisU * Scale

>>> UPDATEVIRIAL
Virial = Virial + ThisW * Scale

>>> UPDATEFORCE
FORCEI = FORCEI * SCALE
FORCE(I,:) = FORCE(I,:) + FORCEI
FORCE(J,:) = FORCE(J,:) - FORCEI
"""

SourceUpdatesNoScale = """
>>> UPDATEENERGY
PEnergy = PEnergy + ThisU
ThisPEnergy = ThisPEnergy + ThisU

>>> UPDATEVIRIAL
Virial = Virial + ThisW

>>> UPDATEFORCE
FORCEI = FORCEI
FORCE(I,:) = FORCE(I,:) + FORCEI
FORCE(J,:) = FORCE(J,:) - FORCEI
"""

SourcePairTemplate = """
>>> argmainlooppair
ARGVAL = DIJ
ARGTYPE = 0
if (ArgVal*ArgVal > CUTSQ) cycle
[ARGGET]

>>> argeval
DIJ = ARGVAL
DIJSQ = DIJ * DIJ
[pairenergy]
[UPDATEENERGY]
if (CALCVIRIAL .or. CALCDWPARAM) then
    [pairvirial]
    [UPDATEVIRIAL]
endif
if (CALCDUPARAM .or. CALCFLUCT) then
    [pairduparam]
endif
if (CALCDWPARAM) then
    [pairdwparam]
endif

>>> mainlooppair
[pairenergy]
[UPDATEENERGY]
if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
    [pairvirial]
    [UPDATEVIRIAL]
endif
if (CALCFORCE) then
    FORCEI = (RIJ * THISW / DIJSQ) * SCALE
    [UPDATEFORCE]
endif
if (CALCDUPARAM .or. CALCFLUCT) then
    [pairduparam]
endif
if (CALCDWPARAM) then
    [pairdwparam]
endif
"""

SourceMain = """
LoopMode = Mode
if (Mode == 0) then
    !zero initial quantities
    PEnergy = 0.d0
    Terms = 0.d0
    if (CalcVirial) then
        Virial = 0.d0
    else
        Virial = 0.d0
    endif
    if (CalcForce) Force = 0.d0
    if (CalcDUParam .or. CalcFluct) then
        DUParam = 0.d0
        DDUParam = 0.d0
    endif
    if (CalcDWParam) then
        DWParam = 0.d0
        DDWParam = 0.d0
    endif
    if (CalcFluct) then
        FluctE0 = 0.
        FluctA0 = 0.
    endif
    Scale = 1.d0
else
    if (CalcDUParam .or. CalcDWParam .or. CalcForce .or. CalcVirial .or. CalcFluct) then
        print *, "Can only use Calcvirial, CalcForce, CalcDUParam, CalcDWParam, CalcFluct for Mode=0."
        stop
    endif
    if (Mode > 0) then
        Scale = 1.d0
    elseif (Mode < 0) then
        Scale = -1.d0
    else
        print *, "Invalid value of Mode in calcenergyforces."
        stop
    endif
endif
"""

SourceSaveEnergyState = """
>>> inputmaps
int Mode = Mode

>>> Main
if (Mode == 0) then
    !save for all atoms
    istart = 0
    istop = NAtom - 1
elseif (Mode == 1) then
    !save for single atom TargetAtom
    istart = TargetAtom
    istop = TargetAtom
elseif (Mode == 2) then
    !save for single molecule TargetMol
    istart = MolRange(TargetMol)
    istop = MolRange(TargetMol+1) - 1
else
    print *, "Illegal value for Mode in SaveEnergyState."
endif
OldPEnergy = PEnergy
OldTerms = Terms
"""

SourceRevertEnergyState = """
>>> inputmaps
int Mode = Mode

>>> Main
if (Mode == 0) then
    !revert for all atoms
    istart = 0
    istop = NAtom - 1
elseif (Mode == 1) then
    !revert for single atom TargetAtom
    istart = TargetAtom
    istop = TargetAtom
elseif (Mode == 2) then
    !revert for single molecule TargetMol
    istart = MolRange(TargetMol)
    istop = MolRange(TargetMol+1) - 1
else
    print *, "Illegal value for Mode in RevertEnergyState."
endif
PEnergy = OldPEnergy
Terms = OldTerms
"""


def ReplaceTokens(Source, FF, P = None, PInd = None, Coef = None):
    """Replaces parameter arrays with appropriate index for potential."""  
    #globals
    for (Name, (Start, Stop)) in FF.Globals.Maps.items():
        Source = fortran.ReplaceVarMaster(Source, Name, "Param",
                                          StartInd = FF.Globals.Ind1,
                                          ShiftInd = (Start,),
                                          Len = Stop - Start)

    #potentials
    if not P is None:
        #now replace parameter names
        for (Name, (Start, Stop)) in P.Param.Maps.items():
            Source = fortran.ReplaceVarMaster(Source, Name, "Param",
                                              StartInd = P.Param.Ind1,
                                              ShiftInd = (Start,),
                                              Len = Stop - Start)
    
    for x in ["U", "W"]:
        
        #FIRST DERIVATIVES
        for (Name, (Start, Stop)) in FF.Globals.Maps.items():
            #check that DParam is not used on RHS of assignment
            DName = "D%s_%s" % (x, Name)
            VName = "D%sParam" % x
            if fortran.HasAssign(Source, DName, DName):
                raise StandardError("%s appears on RHS of assignment in:\n%s" % (DName, Source))
            #replace assignment expressions with augmentations (e.g., += )
            Source = fortran.AugmentAssign(Source, DName, Coef)
            #now search and replace with master DParam array
            Source = fortran.ReplaceVarMaster(Source, DName, VName,
                                              StartInd = FF.Globals.Ind1,
                                              ShiftInd = (Start,),
                                              Len = Stop - Start)
    
        if not P is None:
            for (Name, (Start, Stop)) in P.Param.Maps.items():
                #check that DParam is not used on RHS of assignment
                DName = "D%s_%s" % (x, Name)
                VName = "D%sParam" % x
                if fortran.HasAssign(Source, DName, DName):
                    raise StandardError("%s appears on RHS of assignment in:\n%s" % (DName, Source))
                #replace assignment expressions with augmentations (e.g., += )
                Source = fortran.AugmentAssign(Source, DName, Coef)
                #now search and replace with master DParam array
                Source = fortran.ReplaceVarMaster(Source, DName, VName,
                                                  StartInd = P.Param.Ind1,
                                                  ShiftInd = (Start,),
                                                  Len = Stop - Start)
           
        #SECOND DERIVATIVES
        for (Name1, (Start1, Stop1)) in FF.Globals.Maps.items():
            for (Name2, (Start2, Stop2)) in FF.Globals.Maps.items():
                #check that DDParam is not used on RHS of assignment
                DName = "DD%s_%s_%s" % (x, Name1, Name2)
                VName = "DD%sParam" % x
                if fortran.HasAssign(Source, DName, DName):
                    raise StandardError("%s appears on RHS of assignment in:\n%s" % (DName, Source))
                #replace assignment expressions with augmentations (e.g., += )
                Source = fortran.AugmentAssign(Source, DName, Coef)
                #now search and replace with master DDParam array
                Source = fortran.ReplaceVarMaster(Source, DName, VName,
                                                  Strides = FF.Globals.Strides2,
                                                  StartInd = FF.Globals.Ind2,
                                                  ShiftInd = (Start1, Start2))
    
        if not P is None:
            for (Name1, (Start1, Stop1)) in P.Param.Maps.items():
                for (Name2, (Start2, Stop2)) in P.Param.Maps.items():
                    #check that DDParam is not used on RHS of assignment
                    DName = "DD%s_%s_%s" % (x, Name1, Name2)
                    VName = "DD%sParam" % x
                    if fortran.HasAssign(Source, DName, DName):
                        raise StandardError("%s appears on RHS of assignment in:\n%s" % (DName, Source))
                    #replace assignment expressions with augmentations (e.g., += )
                    Source = fortran.AugmentAssign(Source, DName, Coef)
                    #now search and replace with master DDParam array
                    Source = fortran.ReplaceVarMaster(Source, DName, VName,
                                                      Strides = P.Param.Strides2,
                                                      StartInd = P.Param.Ind2,
                                                      ShiftInd = (Start1, Start2))
                    
            for (Name1, (Start1, Stop1)) in FF.Globals.Maps.items():
                for (Name2, (Start2, Stop2)) in P.Param.Maps.items():
                    #GLOBAL FIRST THEN POTENTIAL
                    #check that DDParam is not used on RHS of assignment
                    DName = "DD%s_%s_%s" % (x, Name1, Name2)
                    VName = "DD%sParam" % x
                    if fortran.HasAssign(Source, DName, DName):
                        raise StandardError("%s appears on RHS of assignment in:\n%s" % (DName, Source))
                    #replace assignment expressions with augmentations (e.g., += )
                    Source = fortran.AugmentAssign(Source, DName, Coef)
                    #now search and replace with master DDParam array
                    Strides = P.Param.Strides2Global
                    Source = fortran.ReplaceVarMaster(Source, DName, VName,
                                                      Strides = Strides,
                                                      StartInd = P.Param.Ind2Global,
                                                      ShiftInd = (Start1, Start2))
                    #POTENTIAL FIRST THEN GLOBAL
                    #enforce global first
                    DName = "DD%s_%s_%s" % (x, Name2, Name1)
                    if fortran.HasToken(Source, DName):
                        raise StandardError("Global %s must appear before %s" % (Name1, Name2))
                                         
    #cutoff
    if not PInd is None:
        Source = fortran.ReplaceToken(Source, "CutSq", "CutSq(%d)" % PInd)   
        Source = fortran.ReplaceToken(Source, "Cut", "Cut(%d)" % PInd)  
        
    #term energy
    if not PInd is None:
        Source = fortran.ReplaceToken(Source, "ThisPEnergy", "Terms(%d)" % PInd)
        
    #module vars, update variables with prefix
    if not P is None:
        for (Var, ModuleVar) in P.LibVars0.items():
            Source = fortran.ReplaceToken(Source, Var, ModuleVar)

    return Source
    

def GetFortCode(FF):
    """Returns defs and fortran code for the potential energy function."""
    #FIRST MAKE THE ENERGY FORCES ROUTINE
    
    #start new fortran code
    Source = SourceMain + fortran.loops.GetSourceVars() + SourceVars
    fc = fortran.FortCode(Source = Source)   
    
    #make force field code objects
    PGroup = []
    PfcTemplate = fortran.FortCode(SourceBase)
    for (PInd,P) in enumerate(FF):
        Pd = {"PName":str(P), "PInd":PInd, "CutSq":"CutSq(%d)" % PInd, "Cut":"Cut(%d)" % PInd}
        #first check code
        fc.CheckFortCode(P.FortCode, ErrOrigin = P.Label)
        #now replace parameters and make new
        Pfc = P.FortCode.Copy()
        if P.NoEnergyUpdateScaling:
            Pfc.Add(SourceUpdatesNoScale)
        else:
            Pfc.Add(SourceUpdates)
        Source = Pfc.RawSource() % Pd
        Source = ReplaceTokens(Source, FF, P = P, PInd = PInd)
        Pfc = fortran.FortCode(Source, Comment = "potential %s" % P.Label)
        Pfc.AddMissing(PfcTemplate)
        Pfc.VarIndex = PInd
        fc.AddFortCodeVars(Pfc)
        PGroup.append((P, Pd, Pfc))
       
    #add the init commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[init]")
        fc.Add(s)
                    
    #make loops                   
    for Prefix in ["preloop", "mainloop"]:
       
        #add the next commands
        for (P,Pd,Pfc) in PGroup:
            s = Pfc.ProcessMain("[before%s]" % Prefix)
            fc.Add(s)
        
        #make loop 
        l = fortran.loops.Loop1or2()
    
        #add commands from each poential
        for (P,Pd,Pfc) in PGroup:
            s = Pfc.ProcessMain("[%sbeforepair]" % Prefix)
            if P.Filter.N == 1: 
                l.AddLoop(SourceSingle = s, Filter = P.Filter, Sys = FF.Sys)
            elif P.Filter.N > 1:
                l.AddLoop(SourceSingle = s, Filter = P.Filter[0], Sys = FF.Sys)
            s = Pfc.ProcessMain("[%spair]" % Prefix)
            if P.Filter.N == 2:
                l.AddLoop(SourcePair = s, Filter = P.Filter, Sys = FF.Sys, PairCutSqStr = Pd["CutSq"])
            s = Pfc.ProcessMain("[%safterpair]" % Prefix)
            if P.Filter.N == 1: 
                l.AddLoop(SourceAfterPair = s, Filter = P.Filter, Sys = FF.Sys)
            elif P.Filter.N > 1:
                l.AddLoop(SourceAfterPair = s, Filter = P.Filter[0], Sys = FF.Sys)
            
        #add loop
        fc.AddFortCode(l) 
        
    fc.Break()
        
    #make local density energy loop; not additive
    l = fortran.loops.Loop1or2(IsAdditive = False)
    for (P,Pd,Pfc) in PGroup:
        if not P.Type == ptypes.LocalDensityPotential: continue
        s = Pfc.ProcessMain("[localdensityloop1]")
        l.AddLoop(SourceSingle = s, Filter = P.Filter[0], Sys = FF.Sys)
    fc.AddFortCode(l)
    
    fc.Break()
    
    #make local density force loop; not additive
    l = fortran.loops.Loop1or2(IsAdditive = False)
    for (P,Pd,Pfc) in PGroup:
        if not P.Type == ptypes.LocalDensityPotential: continue
        s = Pfc.ProcessMain("[localdensityloop2]")
        l.AddLoop(SourcePair = s, Filter = P.Filter, Sys = FF.Sys, PairCutSqStr = Pd["CutSq"])
    fc.AddFortCode(l, IfCondition = "CalcForce .or. CalcVirial .or. CalcDWParam")
       
    #add the next commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[beforeangleloop]")
        fc.Add(s)   
        
    #angle loop
    l = fortran.loops.LoopNIntra(N=3)
    for (P,Pd,Pfc) in PGroup:
        if not P.Type == ptypes.AnglePotential: continue
        s = Pfc.ProcessMain("[angleloop]")
        l.AddLoop(s, Filter = P.Filter, Sys = FF.Sys)
    
    #add angle loop
    fc.AddFortCode(l)
    
    #add the next commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[beforetorsionloop]")
        fc.Add(s)
    
    #torsion loop
    l = fortran.loops.LoopNIntra(N=4)
    for (P,Pd,Pfc) in PGroup:
        if not P.Type == ptypes.TorsionPotential: continue
        s = Pfc.ProcessMain("[torsionloop]")
        l.AddLoop(s, Filter = P.Filter, Sys = FF.Sys)
    
    #add torsion loop
    fc.AddFortCode(l)   

    #add the final commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[final]")
        fc.Add(s)
              
    #NOW MAKE THE SAVE/REVERT ENERY STATE ROUTINES
        
    #start new fortran code
    Source = SourceSaveEnergyState + fortran.loops.GetSourceVars() + SourceVars
    fc1 = fortran.FortCode(Source = Source)   
    Source = SourceRevertEnergyState + fortran.loops.GetSourceVars() + SourceVars
    fc2 = fortran.FortCode(Source = Source) 
    
    fc1.Break()
    fc2.Break()
          
    #add the save/revert energy state commands
    for (P,Pd,Pfc) in PGroup:
        s = Pfc.ProcessMain("[saveenergystate]")
        fc1.Add(s)
        s = Pfc.ProcessMain("[revertenergystate]")
        fc2.Add(s)

    return fc, fc1, fc2
        
        
