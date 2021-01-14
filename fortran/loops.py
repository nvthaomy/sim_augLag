#/usr/bin/env python

### Fortran code generation for loops

import numpy as np

from codeblock import FortCode
from code import IfBracket, HasToken, ReplaceToken, CodeStrip


AtomIndices = ["i", "j", "k", "l"]


SourceVars = r"""
>>> inputmaps

>>> outputmaps

>>> defs
float iBoxL(Dim)
bool DoMinImage
int i
int j
int k
int l
int m
int ForceDoMInd
int ForceDoAtom
int AIDi
int AIDj
int AIDmin
int AIDmax
int AIDij
int SIDi
int SIDj
int SIDmin
int SIDmax
int MIDi
int MIDj
int MIDm
int MIDmin
int MIDmax
int MIndi
int MIndj
int ApInd1
int ApInd2
int ApInd
int BondOrdIndShift
int BondOrdSIDjStart
int BondOrdSIDjStop
int BondOrdij
int istart
int istop
int jstart
int mstart
int mstop
int LoopMode
bool SameMol
bool Bonded
bool TorDefined
bool Applyi
bool Applyj
bool AnyActive
bool MolIsRigidi
bool SameRigidMol
int ListIndi
float Posi(Dim)
float Posj(Dim)
float PosiMinImage(Dim)
float PosjMinImage(Dim)
float Massi
float Massj
float rij(Dim)
float dijsq
float dij
float rjk(Dim)
float djksq
float djk
float rkl(Dim)
float dklsq
float dkl        
float theta
float costheta
float sintheta
float phi
float crossijk(Dim)
float crossjkl(Dim)
float crosssqijk
float crosssqjkl
float pi = 3.141592653589d0
float hpi = 1.570796326795d0
float sqrtpi = 1.772453850905d0
float thetaijk
float costhetaijk
float sinthetaijk
float thetajkl
float costhetajkl
float sinthetajkl
"""
          


SourceFilter1On = """
ForceDoMInd = -1
ForceDoAtom = -1
if (LoopMode == 0) then
    !normal full atom loop
    istart = 0
    istop = NAtom - 1
elseif (LoopMode == 1) then
    !single atom interactions for adding TargetAtom
    istart = TargetAtom
    istop = TargetAtom
    ForceDoAtom = TargetAtom
elseif (LoopMode == -1) then
    !single atom interactions for deleting TargetAtom
    istart = TargetAtom
    istop = TargetAtom
    ForceDoAtom = TargetAtom
elseif (LoopMode == 2) then
    !molecule interactions fo adding TargetMol
    istart = MolRange(TargetMol)
    istop = MolRange(TargetMol+1) - 1
    ForceDoMInd = TargetMol
elseif (LoopMode == -2) then
    !molecule interactions fo deleting TargetMol
    istart = MolRange(TargetMol)
    istop = MolRange(TargetMol+1) - 1
    ForceDoMInd = TargetMol
else
    print *, "Illegal value for LoopMode."
endif
"""

SourceFilter2On = """
if (MolActive(MIndi) < 0 .and. LoopMode == 0) cycle    
"""

SourceFilter3On = """
if (LoopMode == 0) then
    jstart = i+1
else
    jstart = 0
endif
"""

SourceFilter4On = """
if (LoopMode == 2 .or. LoopMode == -2) then
    !!check to see if we need to skip because of double counting
    if (MIndj == MIndi .and. j < i) cycle
endif

if (MolActive(MIndj) < 0 .and. MIndj /= ForceDoMInd .and. j /= ForceDoAtom) cycle
"""

SourceFilter1Off = """
!normal full atom loop
istart = 0
istop = NAtom - 1
"""

SourceFilter2Off = """
if (LoopMode == 0) then
    if (MolActive(MIndi) < 0) cycle
elseif (LoopMode == 1) then
    if (MolActive(MIndi) < 0 .and. TargetAtom /= i) cycle
elseif (LoopMode == -1) then
    if (MolActive(MIndi) < 0 .or. TargetAtom == i) cycle
elseif (LoopMode == 2) then
    if (MolActive(MIndi) < 0 .and. TargetMol /= MIndi) cycle
elseif (LoopMode == -2) then
    if (MolActive(MIndi) < 0 .or. TargetMol == MIndi) cycle
endif
"""

SourceFilter3Off = """
jstart = i+1
"""

SourceFilter4Off = """
if (LoopMode == 0) then
    if (MolActive(MIndj) < 0) cycle
elseif (LoopMode == 1) then
    if (MolActive(MIndj) < 0 .and. TargetAtom /= j) cycle
elseif (LoopMode == -1) then
    if (MolActive(MIndj) < 0 .or. TargetAtom == j) cycle
elseif (LoopMode == 2) then
    if (MolActive(MIndj) < 0 .and. TargetMol /= MIndj) cycle
elseif (LoopMode == -2) then
    if (MolActive(MIndj) < 0 .or. TargetMol == MIndj) cycle
endif
"""
           
Source12 = r"""
DoMinImage = any(BoxL > 0.d0)                                      || DoMinImage
iBoxL = 1.d0 / max(1.d-300, BoxL)                                  || iBoxL

[LOOPFILTER1]                                                      || istart

!loop over i                                                       ||+ [SourceSingle], [SourcePair], [SourceAfterPair]
do i = istart, istop

    MIndi = MInd(i) 
    
    [LOOPFILTER2]                                

    Posi = Pos(i,:)                                                || Posi
    PosiMinImage = Posi                                            || PosiMinImage
    if (DoMinImage) PosiMinImage = Posi - BoxL * dnint(Posi * iBoxL)   || PosiMinImage
    AIDi = AID(i)                                                  || AIDi
    SIDi = SID(i)                                                  || SIDi
    MIDi = MID(i)                                                  || MIDi
    Massi = Mass(i)                                                || Massi
    BondOrdSIDjStart = BondOrdShift(SIDi)                          || BondOrdij
    BondOrdSIDjStop = BondOrdSIDjStart + BondOrdStart(SIDi+1) - BondOrdStart(SIDi) || BondOrdij
    BondOrdIndShift = BondOrdStart(SIDi) - BondOrdShift(SIDi)      || BondOrdij
    MolIsRigidi = MolIsRigid(MIndi)                                || MolIsRigidi
    
    [SourceSingle]
    
    [LOOPFILTER3]                                                  || [SourcePair]
    
    !loop over j                                                   ||+ [SourcePair]
    do j = jstart, NAtom - 1
    
        !check to see if same atom
        if (i==j) cycle

        MIndj = MInd(j) 
        
        [LOOPFILTER4]

        Posj = Pos(j,:)                                            || Posj
        PosjMinImage = Posj                                        || PosjMinImage
        if (DoMinImage) PosjMinImage = Posj - BoxL * dnint(Posj * iBoxL)   || PosjMinImage
        AIDj = AID(j)                                              || AIDj
        SIDj = SID(j)                                              || SIDj
        MIDj = MID(j)                                              || MIDj
        if (SIDi < SIDj) then                                      || AIDmin, AIDmax, SIDmin, SIDmax, MIDmin, MIDmax
            AIDmin = AIDi                                          || AIDmin, AIDmax
            AIDmax = AIDj                                          || AIDmin, AIDmax
            SIDmin = SIDi                                          || SIDmin, SIDmax
            SIDmax = SIDj                                          || SIDmin, SIDmax
            MIDmin = MIDi                                          || MIDmin, MIDmax
            MIDmax = MIDj                                          || MIDmin, MIDmax
        else                                                       || AIDmin, AIDmax, SIDmin, SIDmax, MIDmin, MIDmax
            AIDmin = AIDj                                          || AIDmin, AIDmax
            AIDmax = AIDi                                          || AIDmin, AIDmax
            SIDmin = SIDj                                          || SIDmin, SIDmax
            SIDmax = SIDi                                          || SIDmin, SIDmax
            MIDmin = MIDj                                          || MIDmin, MIDmax
            MIDmax = MIDi                                          || MIDmin, MIDmax
        endif                                                      || AIDmin, AIDmax, SIDmin, SIDmax, MIDmin, MIDmax
        if (AIDi > AIDj) then                                      ||+ AIDij
            AIDij = AIDi * (AIDi + 1) / 2 + AIDj
        else
            AIDij = AIDj * (AIDj + 1) / 2 + AIDi
        endif                                                    
        SameMol = (MIndi == MIndj)                                 || SameMol
        BondOrdij = BondOrdLimit + 1                               || BondOrdij
        if (SameMol) then                                          || BondOrdij 
            if (SIDj >= BondOrdSIDjStart .and. SIDj < BondOrdSIDjStop) then  || BondOrdij
                BondOrdij = BondOrdData(SIDj + BondOrdIndShift)    || BondOrdij
            endif                                                  || BondOrdij
        endif                                                      || BondOrdij
        Bonded = (SameMol .and. BondOrdij==2)                      || Bonded
        Massj = Mass(j)                                            || Massj
        
        rij = Posj - Posi                                          || rij
        if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)      || rij
        dijsq = dot_product(rij, rij)                              || dijsq
        dij = -1.d0                                                || dij
        
        SameRigidMol = (MIndi == MIndj) .and. MolIsRigidi          || SameRigidMol

        [SourcePair]
            
    !end of loop j                                                 || [SourcePair]
    enddo                                                          || [SourcePair]
    
    [SourceAfterPair]

!end of loop i                                                     || [SourceSingle], [SourcePair], [SourceAfterPair]
enddo                                                              || [SourceSingle], [SourcePair], [SourceAfterPair]
"""    


Source3Intra = r"""
DoMinImage = any(BoxL > 0.d0)                                      || DoMinImage
iBoxL = 1.d0 / max(1.d-300, BoxL)                                  || iBoxL

if (LoopMode == 0) then                                                ||+ [Source]
    !normal full atom loop
    mstart = 0
    mstop = NMol - 1
elseif (LoopMode == 1 .or. LoopMode == -1) then
    !single atom interactions
    mstart = MInd(TargetAtom)
    mstop = MInd(TargetAtom)
elseif (LoopMode == 2 .or. LoopMode == -2) then
    mstart = TargetMol
    mstop = TargetMol
endif

!3-atom intramolecular loop                                        ||+ [Source]
do m = mstart, mstop

    if (MolActive(m) < 0 .and. LoopMode == 0) cycle    
    
    MIDm = MolID(m)

    [SourceApply]
    
    do ApInd = ApInd1, ApInd2, 3
        i = %(IndVar)s(ApInd) + MolRange(m)
        j = %(IndVar)s(ApInd+1) + MolRange(m)
        k = %(IndVar)s(ApInd+2) + MolRange(m)
        
        rij = Pos(j,:) - Pos(i,:)                                  || rij
        if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)      || rij
        dij = sqrt(sum(rij*rij))                                   || dij
        rjk = Pos(k,:) - Pos(j,:)                                  || rjk
        if (DoMinImage) rjk = rjk - BoxL * dnint(rjk * iBoxL)      || rjk
        djk = sqrt(sum(rjk*rjk))                                   || djk
        costheta = -dot_product(rij,rjk) / (dij * djk)             || costheta
        costheta = max(min(costheta, 1.d0), -1.d0)                 || costheta
        sintheta = sqrt(max(0.d0, 1.d0 - costheta*costheta))       || sintheta
        theta = acos(costheta)                                     || theta
        
        [Source]                                                   ||
    enddo                                                          || [Source]
enddo                                                              || [Source]
"""


Source4Intra = r"""
DoMinImage = any(BoxL > 0.d0)                                      || DoMinImage
iBoxL = 1.d0 / max(1.d-300, BoxL)                                  || iBoxL

if (LoopMode == 0) then                                                ||+ [Source]
    !normal full atom loop
    mstart = 0
    mstop = NMol - 1
elseif (LoopMode == 1 .or. LoopMode == -1) then
    !single atom interactions
    mstart = MInd(TargetAtom)
    mstop = MInd(TargetAtom)
elseif (LoopMode == 2 .or. LoopMode == -2) then
    mstart = TargetMol
    mstop = TargetMol
endif

!4-atom intramolecular loop                                        ||+ [Source]
do m = mstart, mstop
 
    if (MolActive(m) < 0 .and. LoopMode == 0) cycle

    MIDm = MolID(m)
    
    [SourceApply]
    
    do ApInd = ApInd1, ApInd2, 4
        i = %(IndVar)s(ApInd) + MolRange(m)
        j = %(IndVar)s(ApInd+1) + MolRange(m)
        k = %(IndVar)s(ApInd+2) + MolRange(m)
        l = %(IndVar)s(ApInd+3) + MolRange(m)
        
        rij = Pos(j,:) - Pos(i,:)                                  || rij
        if (DoMinImage) rij = rij - BoxL * dnint(rij * iBoxL)      || rij
        dij = sqrt(sum(rij*rij))                                   || dij
        rjk = Pos(k,:) - Pos(j,:)                                  || rjk
        if (DoMinImage) rjk = rjk - BoxL * dnint(rjk * iBoxL)      || rjk
        djk = sqrt(sum(rjk*rjk))                                   || djk
        rkl = Pos(l,:) - Pos(k,:)                                  || rkl
        if (DoMinImage) rkl = rkl - BoxL * dnint(rkl * iBoxL)      || rkl
        dkl = sqrt(sum(rkl*rkl))                                   || dkl
        costhetaijk = -dot_product(rij,rjk) / (dij * djk)          || costhetaijk
        costhetaijk = max(min(costhetaijk, 1.d0), -1.d0)           || costhetaijk
        sinthetaijk = sqrt(max(0.d0, 1.d0 - costhetaijk*costhetaijk)) || sinthetaijk
        thetaijk = acos(costhetaijk)                               || thetaijk
        costhetajkl = -dot_product(rjk,rkl) / (djk * dkl)          || costhetajkl
        costhetajkl = max(min(costhetajkl, 1.d0), -1.d0)           || costhetajkl
        sinthetajkl = sqrt(max(0.d0, 1.d0 - costhetajkl*costhetajkl)) || sinthetajkl
        thetajkl = acos(costhetajkl)                               || thetajkl
        crossijk = (/ rij(1)*rjk(2) - rij(2)*rjk(1), &             || crossijk
                    & rij(2)*rjk(0) - rij(0)*rjk(2), &             || crossijk
                    & rij(0)*rjk(1) - rij(1)*rjk(0) /)             || crossijk
        crosssqijk = dot_product(crossijk, crossijk)               || crosssqijk
        crossjkl = (/ rjk(1)*rkl(2) - rjk(2)*rkl(1), &             || crossjkl
                    & rjk(2)*rkl(0) - rjk(0)*rkl(2), &             || crossjkl
                    & rjk(0)*rkl(1) - rjk(1)*rkl(0) /)             || crossjkl
        crosssqjkl = dot_product(crossjkl, crossjkl)               || crosssqjkl
        TorDefined = (crosssqijk > 0. .and.  crosssqjkl > 0.)      || TorDefined
        if (TorDefined) then                                       ||+ phi
            phi = dot_product(crossijk, crossjkl) / (sqrt(crosssqijk*crosssqjkl))
            phi = min(1.d0, max(-1.d0, phi))
            phi = sign(acos(phi), dot_product(rij, crossjkl))
        else
            phi = 0.d0
        endif                                                      

        [Source]                                                   ||
    enddo                                                          || [Source]
enddo                                                              || [Source]
"""


def FilterFortran(F, World, Ind = AtomIndices[0]):
    """Returns a fortran boolean expression for a filter in loops."""
    from code import AddParen
    sTypes = []
    if F.HasSiteType():
        SiteTypes = F.Select(World.SiteTypes)
        sTypes.extend(["SID%s==%d" % (Ind, a.SID) for a in SiteTypes])
    if F.HasAtomType():
        AtomTypes = F.Select(World.AtomTypes)
        sTypes.extend(["AID%s==%d" % (Ind, a.AID) for a in AtomTypes])
    sTypes = " .or. ".join(sTypes)
    if not F.MolTypes is None:
        sMolTypes = " .or. ".join(["MID%s==%d" % (Ind, m.MID) for m in F.MolTypes])
    else:
        sMolTypes = ""
    if sTypes and sMolTypes:
        s = "%s .and. %s" % (AddParen(sTypes), AddParen(sMolTypes))
    elif sTypes:
        s = sTypes
    elif sMolTypes:
        s = sMolTypes
    else:
        return ""
    if F.Negate:
        return ".not. %s" % AddParen(s)
    else:
        return s
        
        
def PolyFilterFortran(F, World, Ind = AtomIndices):
    """Returns a fortran boolean expression for loops."""
    from code import AddParen
    import sim.chem as chem
    if F.N > len(Ind):
        raise IndexError("Too few index variables defined.")
    if F.N == 1:
        return FilterFortran(F[0], World, Ind[0])
    #this only works if the length is 2
    if F.N > 2:
        raise ValueError("Cannot run this routine with more than two atom specs in filter.")
    #make a new index list
    Ind = Ind[:F.N]
    #initialize the specification list
    sSpec = []
    #do the atom specifications
    if F.Ordered:
        for (Filter, ThisInd) in zip(F, Ind):
            sSpec.append(FilterFortran(Filter, World, ThisInd))
    else:
        #need to account for forward and reverse
        sSpec1 = []
        sSpec2 = []
        for (Filter, ThisInd) in zip(F, Ind):
            sSpec1.append(FilterFortran(Filter, World, ThisInd))
        for (Filter, ThisInd) in zip(F, Ind[::-1]):
            sSpec2.append(FilterFortran(Filter, World, ThisInd))
        sSpec1 = " .and. ".join([AddParen(s) for s in sSpec1 if len(s)])
        sSpec2 = " .and. ".join([AddParen(s) for s in sSpec2 if len(s)])
        if len(sSpec1):
            sSpec.append(AddParen(sSpec1) + " .or. " + AddParen(sSpec2))
    #check the bonded specifications
    if F.Bonded == True:
        sSpec.append("Bonded")
    elif F.Bonded == False:
        sSpec.append("(.not. Bonded)")
    elif not F.MaxBondOrd is None or not F.MinBondOrd is None:
        if not F.MinBondOrd is None:
            if F.MinBondOrd > chem.BondOrdLimit + 1:
                raise ValueError("MinBondOrd of %d is larger than tabulated limit of %d in chem.BondOrdLimit" % (F.MinBondOrd, chem.BondOrdLimit))
            sSpec.append("BondOrdij>=%d" % F.MinBondOrd)
        if not F.MaxBondOrd is None:
            if F.MaxBondOrd > chem.BondOrdLimit:
                raise ValueError("MaxBondOrd of %d is larger than tabulated limit of %d in chem.BondOrdLimit" % (F.MaxBondOrd, chem.BondOrdLimit))
            sSpec.append("BondOrdij<=%d" % F.MaxBondOrd)
    elif F.Intra == True:
        sSpec.append("SameMol")
    elif F.Intra == False:
        sSpec.append("(.not. SameMol)")
    #clean up the specifications
    sSpec = [s for s in sSpec if len(s)]        
    if len(sSpec) == 0:
        return ""
    elif len(sSpec) == 1:
        s = sSpec[0]
    else:
        sSpec = ["%s" % AddParen(s) for s in sSpec]
        s = " .and. ".join(sSpec)
    return s
    
    
def FilterListFortran(F, World, Ind = AtomIndices):
    """Returns a fortran boolean expression for loops."""
    from code import AddParen
    if F.N > len(Ind):
        raise IndexError("Too few index variables defined.")
    sSpec = [FilterFortran(Filter, World, Ind) for Filter in F]
    #clean up the specifications
    sSpec = [s for s in sSpec if len(s)]        
    if len(sSpec) == 0:
        return ""
    elif len(sSpec) == 1:
        return sSpec[0]
    else:
        sSpec = ["%s" % AddParen(s) for s in sSpec]
        return " .or. ".join(sSpec)
        
        
def GenericFilterFortran(F, World, Ind = AtomIndices):
    import sim.atomselect as atomselect
    if isinstance(F, atomselect.Filter):
        return FilterFortran(F, World, Ind[0])
    elif isinstance(F, atomselect.PolyFilter):
        return PolyFilterFortran(F, World, Ind)
    elif isinstance(F, atomselect.FilterList):
        return FilterListFortran(F, World, Ind)
    else:
        raise TypeError("Don't recognize filter type of %s." % type(F))



def GetSourceVars(SysStr = None):
    """Returns variables."""
    if not SysStr is None:
        SysStr = SysStr.strip()
        return ReplaceToken(SourceVars, "Sys", SysStr)
    else:
        return SourceVars
    


class Loop1or2(FortCode):
    
    def __init__(self, SysStr = None, ApplyDefPrefix = None, IsAdditive = True, IfCondition = None):
        """Initializes a new single or atom pair loop."""
        FortCode.__init__(self, Blocks = ["SourceSingle", "SourcePair", "SourceAfterPair"])
        FortCode.Add(self, GetSourceVars(SysStr))
        #add base code
        self.Blocks["Main"].Add(Source12)
        if ApplyDefPrefix is None:
            self.ApplyDefPrefix = "Apply2"
        else:
            self.ApplyDefPrefix = ApplyDefPrefix
        self.NApplyDef = 0   
        if IsAdditive:
            FortCode.Add(self, SourceFilter1On, Block = "LOOPFILTER1")
            FortCode.Add(self, SourceFilter2On, Block = "LOOPFILTER2")
            FortCode.Add(self, SourceFilter3On, Block = "LOOPFILTER3")
            FortCode.Add(self, SourceFilter4On, Block = "LOOPFILTER4")
        else:
            FortCode.Add(self, SourceFilter1Off, Block = "LOOPFILTER1")
            FortCode.Add(self, SourceFilter2Off, Block = "LOOPFILTER2")    
            FortCode.Add(self, SourceFilter3Off, Block = "LOOPFILTER3") 
            FortCode.Add(self, SourceFilter4Off, Block = "LOOPFILTER4") 
        self.IfCondition = IfCondition
        
    def __NewApplyDef(self, Filter, Sys):
        """Creates a new constant variable with true/false values for each AID pair type. 
Works when Filter.AIDPairs is present.  Returns IndVar."""
        from sim.chem import NPairID, PairID
        if Filter.AIDPairs is None:
            return None
        NAID = Sys.World.NAID
        if not NAID == len(Filter.AIDPairs):
            raise ValueError("AIDPairs is not the same shape as Sys.World.NAID")
        n = NPairID(NAID)
        a = np.zeros(n, dtype=bool)
        for i in range(NAID):
            for j in range (i, NAID):
                a[PairID(i,j)] = Filter.AIDPairs[i,j]                   
        IndVar = self.ApplyDefPrefix + "_%d" % self.NApplyDef
        self.NApplyDef += 1
        Val = "(/ " + ",".join(["%d" % x for x in a]) + " /)"
        Def = "bool %s(%d) = %s" % (IndVar, n, Val)
        self.Defs.append(Def)
        return IndVar        
        
    def AddLoop(self, SourceSingle = "", SourcePair = "", SourceAfterPair = "", 
                Filter = None, Sys = None, 
                PairCutSqStr = None, IfCondition = None,
                Check = False):
        """Adds source to the code block."""
        #check filter
        if not Filter is None and Sys is None:
            raise ValueError("Must specify Sys along with Filter")
        #check for single code
        SourceSingle = CodeStrip(SourceSingle)
        if len(SourceSingle):
            if not Filter is None:
                if not Filter.N == 1:
                    raise TypeError("Filter is not for single atoms but found SourceSingle.")
                Condition = GenericFilterFortran(Filter, Sys.World)
                SourceSingle = IfBracket(SourceSingle, Condition)
            if not IfCondition is None:
                SourceSingle = IfBracket(SourceSingle, IfCondition)
            self.SetCurrent("SourceSingle")
            FortCode.Add(self, SourceSingle, Check=Check)
            self.Break()
        #check for pair code
        SourcePair = CodeStrip(SourcePair)
        if len(SourcePair):
            if not Filter is None:
                if not Filter.N == 2:
                    raise TypeError("Filter is not for atom pairs but found SourcePair.")
                #check specialized pair filter
                IndVar = self.__NewApplyDef(Filter, Sys)
                if Filter.Ordered:
                    Conditioni = GenericFilterFortran(Filter, Sys.World)
                    Conditionj = GenericFilterFortran(Filter.Reverse(), Sys.World)
                    Source1 = ""
                    if len(Conditioni):
                        Source1 += "Applyi = (%s)\n" % Conditioni
                    else:
                        Source1 += "Applyi = .true.\n"
                    if len(Conditionj):
                        Source1 += "Applyj = (%s)\n" % Conditionj
                    else:
                        Source1 += "Applyj = .true.\n"   
                    if IndVar is None:
                        SourcePair = IfBracket(SourcePair, "Applyi .or. Applyj")
                    else:
                        SourcePair = IfBracket(SourcePair, "%s(AIDIJ) .and. (Applyi .or. Applyj)" % IndVar)
                    SourcePair = Source1 + SourcePair
                else:
                    Condition = GenericFilterFortran(Filter, Sys.World)
                    if IndVar is None:
                        SourcePair = IfBracket(SourcePair, Condition)
                    else:
                        if len(Condition):
                            SourcePair = IfBracket(SourcePair, "%s(AIDIJ) .and. (%s)" % (IndVar, Condition))
                        else:
                            SourcePair = IfBracket(SourcePair, "%s(AIDIJ)" % (IndVar))
            if HasToken(SourcePair, "dij"):
                SourcePair = "if (dij < 0.d0) dij = sqrt(dijsq)\n" + SourcePair
            if not PairCutSqStr is None:
                SourcePair = IfBracket(SourcePair, "dijsq < %s" % PairCutSqStr)
            if not IfCondition is None:
                SourcePair = IfBracket(SourcePair, IfCondition)
            self.SetCurrent("SourcePair")
            FortCode.Add(self, SourcePair, Check=Check)
            self.Break() 
        #check for after pair code
        SourceAfterPair = CodeStrip(SourceAfterPair)
        if len(SourceAfterPair):
            if Check: 
                self.CheckSource(SourceAfterPair)
            if not Filter is None:
                if not Filter.N == 1:
                    raise TypeError("Filter is not for single atoms but found SourceAfterPair.")
                Condition = GenericFilterFortran(Filter, Sys.World)
                SourceAfterPair = IfBracket(SourceAfterPair, Condition)
            if not IfCondition is None:
                SourceAfterPair = IfBracket(SourceAfterPair, IfCondition)
            self.SetCurrent("SourceAfterPair")
            FortCode.Add(self, SourceAfterPair, Check=Check)
            self.Break()
                
        
class LoopNIntra(FortCode):
    def __init__(self, N, SysStr = "Sys", Indent = 0, ApplyDefPrefix = None):
        """Initializes a new multi-atom INTRAmolecular loop.
N gives the number of atoms in each returned group."""
        if not (N==3 or N==4):
            raise ValueError("Can only do N=3 or N=4.")
        FortCode.__init__(self)
        self.N = N
        if ApplyDefPrefix is None:
            self.ApplyDefPrefix = "Apply%d" % N
        else:
            self.ApplyDefPrefix = ApplyDefPrefix
        self.NApplyDef = 0     
        #add maps
        FortCode.Add(self, GetSourceVars(SysStr))
        
    def __NewApplyDef(self, Filter, Sys):
        """Creates a new constant variable containing list of AInd
for different molecules. Returns IndVar, IndRange.
IndRange: list of (MID,Start,Stop) index in Def"""
        #find all combinations of atom sitetypes in the world that 
        #work with this filter and organize by aind's within molecule number
        Groups = Filter.Select(Sys.World.SiteTypes)
        d = {}
        for Group in Groups:
            AIndList = [a.AInd for a in Group]
            MID = Group[0].MID
            if not MID in d:
                d[MID] = []
            d[MID].extend(AIndList)
        l = d.items()
        l.sort()
        #make a variable containing all aind's
        AllInd = []
        IndRange = []
        for (MID, AIndList) in l:
            ApInd1 = len(AllInd)
            ApInd2 = ApInd1 + len(AIndList)
            IndRange.append((MID, ApInd1, ApInd2))
            AllInd.extend(AIndList)
        if len(AllInd):
            #make a constant variable name for holding the atom indices
            IndVar = self.ApplyDefPrefix + "_%d" % self.NApplyDef
            self.NApplyDef += 1
            Val = "(/ " + ",".join(["%d" % x for x in AllInd]) + " /)"
            Def = "int %s(%d) = %s" % (IndVar, len(AllInd), Val)
            self.Defs.append(Def)
            return IndVar, IndRange
        else:
            return None, None

    def AddLoop(self, Source = "", Filter = None, Sys = None, 
                IfCondition = None, Check = False):
        """Adds source to the code block."""
        #add code
        Source = CodeStrip(Source)
        if len(Source) == 0:
            return
        #check conditions on Filter
        if Filter is None or Sys is None:
            raise ValueError("Must specify Filter and Sys.")
        elif not Filter.N == self.N:
            raise ValueError("Filter must specify %d atoms." % self.N)
        elif Filter.Bonded == False and Filter.Intra == False:
            raise ValueError("Filter must be intramolecular.")
        #make a new constant apply def
        IndVar, IndRange = self.__NewApplyDef(Filter, Sys)
        if IndVar == None and IndRange == None:
            print "WARNING: Filter %s does not evaluate to any atom types." % str(Filter)
            return
            #raise ValueError("Filter %s does not evaluate to any atom types." % str(Filter))
        #now make the apply code
        fc = FortCode(Blocks = ["SourceApply", "Source"])
        d = {"IndVar" : IndVar, "N" : self.N}
        if self.N == 3:
            fc.Add(Source3Intra % d, Block = "Main")
        elif self.N == 4:
            fc.Add(Source4Intra % d, Block = "Main")
        #select range of indices for this MID
        fc.SetCurrent("SourceApply")
        fc.Add("ApInd1 = -1\nApInd2 = -2", Dep = "[Source]")
        First = True
        for (MID, ApInd1, ApInd2) in IndRange:
            if First:
                fc.Add("if (MIDm==%d) then" % MID, Dep = "[Source]")
                First = False
            else:
                fc.Add("elseif (MIDm==%d) then" % MID, Dep = "[Source]")
            fc.Add("ApInd1 = %d\nApInd2 = %d\n" % (ApInd1, ApInd2-1),
                   Dep = "[Source]", Indent = 4)
        fc.Add("endif", Dep = "[Source]")
        fc.Add("if (ApInd1 < 0) cycle", Dep = "[Source]")
        #add the actual source
        fc.Add(Source, Check=Check, Block = "Source")
        #check if condition
        if not IfCondition is None:
            fc.IfCondition = IfCondition
        #now add all of the code
        FortCode.AddFortCode(self, fc)
        self.Break()
      
        
