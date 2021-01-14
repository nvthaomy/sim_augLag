#/usr/bin/env python

import numpy as np
import base
import sim.chem as chem
from math import erfc

pi = np.pi
Const1 = 2./np.sqrt(pi)


class Ewald(base.PotentialClass):
    """Ewald long-range electrostatic interactions."""

    Names = ["ewald"]

    Source = """
>>> externals
erfc

>>> defs
float Chargei
float Chargej
float idist
float val1
float val2
float val3
float val4
float val5
float erfcterm
float Temp1
float Temp2
float DipoleSum
float OldMinL
Float OldBoxL(Dim)
float Temp3(Dim)
float Temp4(Dim)
float EWScale
int KMaxSq
int KVecSq
int TotK
int i
int j
int k
complex eikri
complex eikrj
  
>>> init
iBoxL = 1.d0 / BoxL
[ewaldcheck]
if (Mode == 0) then
    !normal full atom loop
    istart = 0
    istop = NAtom - 1
    EWScale = 1.
elseif (Mode == 1) then
    !single atom interactions for adding TargetAtom
    istart = TargetAtom
    istop = TargetAtom
    EWScale = 1.
elseif (Mode == -1) then
    !single atom interactions for deleting TargetAtom
    istart = TargetAtom
    istop = TargetAtom
    EWScale = -1.
elseif (Mode == 2) then
    !molecule interactions for adding TargetMol
    istart = MolRange(TargetMol)
    istop = MolRange(TargetMol+1) - 1
    EWScale = 1.
elseif (Mode == -2) then
    !molecule interactions for deleting TargetMol
    istart = MolRange(TargetMol)
    istop = MolRange(TargetMol+1) - 1
    EWScale = -1.
endif

[ewaldupdatevars]
[ewaldrecip]

>>> ewaldcheck
if (any(BOXL /= EWBoxL)) then
    OldMinL = EWMinL
    OldBoxL = EWBoxL
    EWBoxL = BOXL
    EWMinL = minval(EWBoxL)
    Alpha = AlphaL/EWMinL
    Term1 = 2. * Alpha / sqrtpi
    Term2 = -Alpha * Alpha
    Temp1 = EWMinL/OldMinL
    !update the truncation options, and the cutoff distance
    call erfc(Alpha * Cut, erfcterm)
    ShiftTerm = erfcterm / Cut
    !check to see if we need to recalculate kvectors, 
    !or can just scale them for an isotropic volume change
    Temp3 = EWBoxL / OldBoxL
    if (abs(Temp3(0)-Temp3(1)) < 1.d-5 .and. abs(Temp3(1)-Temp3(2)) < 1.e-5 .and. all(OldBoxL > 0.)) then
        KSq = KSq / (Temp1 * Temp1)
        KUTerm = KUTerm * Temp1 * Temp1
        KFTerm = KFTerm / (Temp1 * Temp1)
    else
        [ewaldsetupkvec]
    endif
endif

>>> ewaldsetupkvec
TotK = 0
K_x = 0.
K_y = 0.
K_z = 0.
KSq = 0.
KUTerm = 0.
KFTerm = 0.
KMaxSq = KMax * KMax
do i = 0, KMax
    do j = -KMax, KMax
        do k = -KMax, KMax
            KVecSq = i*i+j*j+k*k
            if ( i==0 ) then
                Temp1 = 1.
            else
                Temp1 = 2.
            endif
            if (KVecSq/=0 .and. KVecSq<=KMaxSq) then
                K_x(TotK) = i
                K_y(TotK) = j
                K_z(TotK) = k
                !calculate k^2
                KSq(TotK) = 4.*Pi*Pi*sum( (real((/i, j, k/)) * iBoxL)**2 )
                !calculate exp(-k^2/4a^2)/k^2
                KUTerm(TotK) = Temp1 * exp(-KSq(TotK)/(4. * Alpha * Alpha))/KSq(TotK)    
                !calculate (4pi/V)*kxyz*exp(-k^2/4a^2)/k^2 for forces
                KFTerm(:,TotK) = 8.*Pi*Pi * Coef * KUTerm(TotK) * (product(iBoxL) * iBoxL)
                KFTerm(0,TotK) = KFTerm(0,TotK) * real(K_x(TotK))
                KFTerm(1,TotK) = KFTerm(1,TotK) * real(K_y(TotK))
                KFTerm(2,TotK) = KFTerm(2,TotK) * real(K_z(TotK))
                TotK = TotK + 1
            endif
        enddo
    enddo
enddo
if (.not. TotK == NKvec) then
    print *, "Total number of K vectors is unexpected in Ewald sum."
    stop
endif

>>> ewaldupdatevars
!calculate e^(ik*r) for kx = 0,1 and ky/kz = -1,0,1
Temp4 = 2.*Pi * iBoxL
do i = iStart, iStop
    !check if mol active
    if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
    Temp3 = -Pos(i,:) * Temp4
    eik_y(KMax-1,i) = dcmplx(cos(Temp3(1)), sin(Temp3(1)))
    eik_z(KMax-1,i) = dcmplx(cos(Temp3(2)), sin(Temp3(2)))
    eik_x(0,i) = dcmplx( 1., 0. )
    eik_y(KMax,i) = dcmplx( 1., 0. )
    eik_z(KMax,i) = dcmplx( 1., 0. )
    Temp3 = Pos(i,:) * Temp4
    eik_x(1,i) = dcmplx(cos(Temp3(0)), sin(Temp3(0)))
    eik_y(KMax+1,i) = dcmplx(cos(Temp3(1)), sin(Temp3(1)))
    eik_z(KMax+1,i) = dcmplx(cos(Temp3(2)), sin(Temp3(2)))
enddo
!calculate remaining e^(ik*r) by recursion
do i = iStart, iStop
    if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
    do j = 2, KMax
        !check if mol active
        eik_x(j,i) = eik_x(1,i) * eik_x(j-1,i)
        eik_y(KMax+j,i) = eik_y(KMax+1,i) * eik_y(KMax+j-1,i)
        eik_z(KMax+j,i) = eik_z(KMax+1,i) * eik_z(KMax+j-1,i)
        eik_y(KMax-j,i) = dconjg(eik_y(KMax+j,i))
        eik_z(KMax-j,i) = dconjg(eik_z(KMax+j,i))
    enddo
enddo
!calculate the structure factor p(k) = sum( e^(ik*r) * q )
if (Mode == 0) then
    Struc = ( 0., 0. )
endif
do i = iStart, iStop
    !check if mol active
    if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
    Chargei = Charge(AID(i))
    do j = 0, NKVec - 1
        eikri = eik_x(K_x(j),i) * eik_y(K_y(j)+KMax,i) * eik_z(K_z(j)+KMax,i)
        Struc(j) = Struc(j) + eikri * Chargei * EWScale
    enddo
enddo

>>> ewaldrecip
Temp2 = 0.5/(Alpha*Alpha)
!val1 = reciprocal energy, val2 = recip virial, val3 = self energy
val1 = 0.
val2 = 0.
do i = 0, NKVec-1
    Temp1 = KUTerm(i) * real(Struc(i) * dconjg(Struc(i)))
    val1 = val1 + Temp1
    val2 = val2 - Temp1 * ( 1. - Temp2 * KSq(i) )
enddo
val1 = val1 * 2. * Pi * Coef * product(iBoxL)
val2 = val2 * 2. * Pi * Coef * product(iBoxL)
!potential self part
val3 = 0.
do i = iStart, iStop
    if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
    Chargei = Charge(AID(i))
    val3 = val3 + Chargei**2
enddo
val3 = -Alpha / sqrtpi * Coef * val3 * EWScale 
if (Mode == 0) then
    !a full update
    ThisU = val1 + val3
    ThisW = val2
    RPEnergy = val1
    RVirial = val2
    SPenergy = val3
else
    !a partial update (ThisU, ThisW represent delta quantities) 
    ThisU = (val1 - RPEnergy) + val3 
    ThisW = val2 - RVirial
    RPEnergy = val1
    RVirial = val2
    SPenergy = SPenergy + val3 
endif
PENERGY = PENERGY + ThisU
THISPENERGY = THISPENERGY + ThisU
if (CALCVIRIAL) VIRIAL = VIRIAL + ThisW
if (CALCFORCE) then
    !forces reciprocal part
    do i = iStart, iStop
        if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
        Temp3 = 0.
        do j = 0, NKVec-1
            eikri = eik_x(K_x(j), i) * eik_y(K_y(j)+KMax, i) * eik_z(K_z(j)+KMax, i)
            Temp3 = Temp3 - KFTerm(:,j) * dimag(dconjg(eikri) * Struc(j))
        enddo
        Temp3 = Temp3 * Charge(AID(i))
        Force(i,:) = Force(i,:) + Temp3 * Scale
    enddo
endif
if (CALCDUPARAM .or. CALCDWPARAM) then
    val1 = 2. * Pi * product(iBoxL) 
    val2 = Coef * val1 
    val3 = -2. * Alpha / sqrtpi
    Temp2 = 0.5/(Alpha*Alpha)
    DU_Coef(0) = ThisU / Coef
    if (CALCDWPARAM) then
        DW_Coef(0) = ThisW / Coef
    endif
    do i = iStart, iStop
        if (MolActive(MInd(i)) < 0 .and. Mode == 0) cycle
        AIDI = AID(i)    
        !reciprocal contribution
        do k = 0, NKVec-1
            val5 = -(1. - Temp2 * KSq(k))
            eikri = eik_x(K_x(k),i) * eik_y(K_y(k)+KMax,i) * eik_z(K_z(k)+KMax,i)
            eikri = dconjg(eikri)
            val4 = 2. * KUTerm(k) * real(Struc(k) * eikri)
            DU_Charge(AIDI) = val2 * val4
            DDU_Charge_Coef(AIDI,0) = val1 * val4
            if (CALCDWPARAM) then
                DW_Charge(AIDI) = val2 * val4 * val5
                DDW_Charge_Coef(AIDI,0) = val1 * val4 * val5
            endif
            do j = i, iStop
                AIDJ = AID(j)
                if (i == j) then
                    Temp1 = 2.
                elseif (AIDI == AIDJ) then
                    Temp1 = 4.
                else
                    Temp1 = 2.
                endif
                eikrj = eik_x(K_x(k),j) * eik_y(K_y(k)+KMax,j) * eik_z(K_z(k)+KMax,j)
                val4 = Temp1 * val2 * KUTerm(k) * real(eikrj * eikri)
                DDU_Charge_Charge(AIDI,AIDJ) = val4
                if (CALCDWPARAM) then
                    DDW_Charge_Charge(AIDI,AIDJ) = val4 * val5
                endif  
            enddo
        enddo
        !self contribution
        DU_Charge(AIDI) = Coef * val3 * Charge(AIDI)
        DDU_Charge_Charge(AIDI, AIDI) = Coef * val3 
        DDU_Charge_Coef(AIDI,0) = val3 * Charge(AIDI)
    enddo
endif

>>> ewaldreal
call erfc(Alpha * DIJ, erfcterm)
if (BondOrdij > 0 .and. BondOrdij <= ExcludeBondOrd) then
    !exclude this interaction so subtract off 1/r
    erfcterm = erfcterm - 1.0
endif
idist = 1.d0 / DIJ
erfcterm = erfcterm * idist
val1 = Coef(0) * Chargei * Chargej 
val2 = erfcterm - ShiftTerm
THISU = val1 * val2 
[UPDATEENERGY]
if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
    val3 = -(erfcterm + Term1 * exp(Term2 * DIJSQ))
    THISW = val1 * val3 
    [UPDATEVIRIAL]
    if (CALCFORCE) then
        FORCEI = (val1 * val3 * idist*idist) * RIJ 
        FORCE(I,:) = FORCE(I,:) + FORCEI
        FORCE(J,:) = FORCE(J,:) - FORCEI
    endif
endif
if (CALCDUPARAM) then
    DU_Coef(0) = Chargei * Chargej * val2
    if (AIDI==AIDJ) then
        DU_Charge(AIDI) = 2.d0 * Coef(0) * Chargej * val2
        DDU_Charge_Coef(AIDI,0) = 2.d0 * Chargej * val2
        DDU_Charge_Charge(AIDI,AIDJ) = 2.d0 * Coef(0) * val2
    else
        DU_Charge(AIDI) = Coef(0) * Chargej * val2
        DU_Charge(AIDJ) = Coef(0) * Chargei * val2
        DDU_Charge_Coef(AIDI,0) = Chargej * val2
        DDU_Charge_Coef(AIDJ,0) = Chargei * val2
        DDU_Charge_Charge(AIDI,AIDJ) = Coef(0) * val2
    endif
endif
if (CALCDWPARAM) then
    DW_Coef(0) = Chargei * Chargej * val3
    if (AIDI==AIDJ) then
        DW_Charge(AIDI) = 2.d0 * Coef(0) * Chargej * val3
        DDW_Charge_Coef(AIDI,0) = 2.d0 * Chargej * val3
        DDW_Charge_Charge(AIDI,AIDJ) = 2.d0 * Coef(0) * val3
    else
        DW_Charge(AIDI) = Coef(0) * Chargej * val3
        DW_Charge(AIDJ) = Coef(0) * Chargei * val3
        DDW_Charge_Coef(AIDI,0) = Chargej * val3
        DDW_Charge_Coef(AIDJ,0) = Chargei * val3
        DDW_Charge_Charge(AIDI,AIDJ) = Coef(0) * val3
    endif
endif

>>> mainloopbeforepair
Chargei = Charge(AIDI)

>>> mainlooppair
!pair interactions for potential %(PName)s
Chargej = Charge(AIDJ)
[ewaldreal]

>>> saveenergystate
OldRPEnergy = RPEnergy
OldRVirial = RVirial
OldSPEnergy = SPEnergy
OldStruc(istart:istop) = Struc(istart:istop)
oldeik_x(:,istart:istop) = eik_x(:,istart:istop)
oldeik_y(:,istart:istop) = eik_y(:,istart:istop)
oldeik_z(:,istart:istop) = eik_z(:,istart:istop)

>>> revertenergystate
RPEnergy = OldRPEnergy
RVirial = OldRVirial
SPEnergy = OldSPEnergy
Struc(istart:istop) = OldStruc(istart:istop)
eik_x(:,istart:istop) = oldeik_x(:,istart:istop)
eik_y(:,istart:istop) = oldeik_y(:,istart:istop)
eik_z(:,istart:istop) = oldeik_z(:,istart:istop)

>>> argmainlooppair
ARGVAL = DIJ
ARGTYPE = AIDIJ
[ARGGET]

>>> argevalconst
AIDIJ = ArgType
call GetijFromPairInd(AIDIJ, AIDi, AIDj)
Chargei = Charge(AIDI)
Chargej = Charge(AIDJ)

>>> argeval
print *, "Evaluation of energy from pair distance distributions does not work for Ewald sum."
stop
"""

    TestArgs = {"Cut":2.5, "Shift":True, "ExcludeBondOrd":2.}
    Type = base.ptypes.PairPotential
    UsesCharge = True
    UsesATypes = True
    NoHistEval = True

    def __init__(self, Sys, Label = "", Cut = None, Filter = None, 
                 ExcludeBondOrd = 0,
                 Fixed = False, FixedCoef = True, EwaldNAtom = None, 
                 AlphaL = 5.6, KMax = 5, Coef = None, Shift = False):
        """Initializes a screened Coulomb potential.
           MN: EwaldNAtom: number of atoms to include in ewald calculations, for testing speed up if only including charged atoms"""
        if not Filter is None:
            print "WARNING: Ewald does not take a Filter. Use ExcludeBondOrd instead."
        if ExcludeBondOrd > 0:
            print "WARNING: Make sure Ewald cutoff is enough to include atoms with bond orders <= %d" % ExcludeBondOrd
        from sim.atomselect import Pairs
        NType, TypeAIDs, TypeLabels = Sys.World.GetPairTypes()
        base.PotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                     Filter = Pairs,
                                     NType = NType, TypeLabels = TypeLabels,
                                     Source = Ewald.Source)
        #set default EwaldNAtom to number of atoms in the system (MN)
        if EwaldNAtom == None:
            self.EwaldNAtom = self.Sys.NAtom
        else:
            self.EwaldNAtom = EwaldNAtom
        #check dimensions
        if any(Sys.BoxL <= 0.) or len(Sys.BoxL) != 3:
            raise ValueError("Ewald can only be used with 3D periodic systems.")
            
        #add module vars
        self.ModuleVars += ["ExcludeBondOrd", "EWBoxL", "EWMinL", "NKVec", "KMax", "AlphaL", "Alpha", 
                            "KSq", "KUTerm", "KFTerm", "K_x", "K_y", "K_z",
                            "ShiftTerm", "Term1", "Term2", "Struc", "OldStruc", 
                            "RPEnergy", "RVirial", "SPEnergy", 
                            "OldRPEnergy", "OldRVirial", "OldSPEnergy",
                            "eik_x", "eik_y", "eik_z", "oldeik_x", "oldeik_y", "oldeik_z"]            

        self.Shift = Shift
        #must fix coef if charges vary
        self.FixedCoef = FixedCoef

        if Coef is None:
            Coef = 1. / Sys.Units.FPE
            
        self.Param.Add("Coef", 1, Value = Coef, Min = 0., Fixed = Fixed or FixedCoef,
                       Scale = 1. / self.Sys.Units.FPE)
                 
        self.KMax = KMax
        self.AlphaL = AlphaL
        self.ExcludeBondOrd = ExcludeBondOrd
                
        #find the maximum number of k vectors
        self.NKVec = 0
        KMaxSq = self.KMax**2
        for i in range(0, KMax+1):
            for j in range(-KMax, KMax+1):
                for k in range(-KMax, KMax+1):
                    KVecSq = i*i+j*j+k*k
                    if KVecSq!=0 and KVecSq<=KMaxSq: 
                        self.NKVec = self.NKVec + 1
                        
        #initialize k vector arrays
        n = self.NKVec
        self.Struc = np.zeros(n, dtype=complex)
        self.OldStruc = np.zeros(n, dtype=complex)
        self.K_x = np.zeros(n, dtype=int)
        self.K_y = np.zeros(n, dtype=int)
        self.K_z = np.zeros(n, dtype=int)
        self.KSq = np.zeros(n, dtype=float)
        self.KUTerm = np.zeros(n, dtype=float)
        self.KFTerm = np.zeros((self.Sys.Dim, n), dtype=float)
        self.EWBoxL = np.zeros_like(Sys.BoxL)
        self.EWMinL = 0.
        self.Alpha = 0.
        self.ShiftTerm = 0.
        self.Term1 = 0.
        self.Term2 = 0.
        self.RPEnergy = 0.
        self.SPEnergy = 0.
        self.RVirial = 0.
        self.OldRPEnergy = 0.
        self.OldSPEnergy = 0.
        self.OldRVirial = 0.
              
        self.Update()
        self.SetTypeInd()
          

    def SetParam(self, Coef = None):
        """Sets parameters."""
        if not Coef is None:
            if self.FixedCoef:
                print "WARNING: Coefficient is fixed.  Not setting."
            else:
                self.Coef = Coef
                self.Update()
        
    def GetTypeInd(self, Type1, Type2):
        """Returns the type number for this potential."""
        return chem.PairID(Type1.AID, Type2.AID)
        
    def SetTypeInd(self, TypeInd = None):
        """Sets the type to be used in Val and DVal, or None for generic."""
        if TypeInd is None:
            self.CurrentChargeCoef = self.Coef * self.Sys.Units.ECharge**2
        else:
            AID1, AID2 = chem.GetPair(TypeInd)
            Charge1 = self.Sys.World.AtomTypes[AID1].Charge 
            Charge2 = self.Sys.World.AtomTypes[AID2].Charge 
            self.CurrentChargeCoef = self.Coef * Charge1 * Charge2      

    def Val(self, x):
        """Returns the value of the potential."""
        if x >= self.Cut:
            return 0.
        return self.CurrentChargeCoef * (1./x - 1./self.Cut)
            
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        return -self.CurrentChargeCoef / x**2        

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.Coef.SoftMin = 0.0001 / self.Sys.Units.FPE
        if not MaxFracChange is None:
            self.Coef.MaxChange = MaxFracChange / self.Sys.Units.FPE       
            
    def PreLoad(self):
        #make the local density arrays
        NAtom = self.EwaldNAtom
        n = self.KMax
        self.eik_x = np.zeros((n+1, NAtom), dtype = complex)
        self.oldeik_x = np.zeros((n+1, NAtom), dtype = complex)       
        self.eik_y = np.zeros((2*n+1, NAtom), dtype = complex)
        self.oldeik_y = np.zeros((2*n+1, NAtom), dtype = complex)        
        self.eik_z = np.zeros((2*n+1, NAtom), dtype = complex)
        self.oldeik_z = np.zeros((2*n+1, NAtom), dtype = complex)        
        base.PotentialClass.PreLoad(self)                  

    def Update(self):
        """Updates the potential."""
        pass #all updates done in fortran

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        Coef = None
        if not self.FixedCoef: Coef = 1. / self.Sys.Units.FPE
        self.SetParam(Coef)
        


 
