#/usr/bin/env python

### Basic velocity verlet integration in SIM suite.


#TODO:
#-enable RATTLE for fixed atoms (i.e., with AtomList[i]=False)
#-currently Nose Hoover does not work for fixed atoms either


from sim.fortran import FortCode


#thermostat codes
Thermostats = {"ThermostatNone" : 0,
               "ThermostatAndersen" : 1,
               "ThermostatAndersenMassive" : 1,
               "ThermostatAndersenCollision" : 2,
               "ThermostatNoseHoover" : 3,
               "ThermostatLangevin" : 4,}

Barostats = {"BarostatNone" : 0,
             "BarostatMonteCarlo" : 1,}
             
#whether or not the momentum is a constant of motion; used for counting DOFs
ThermostatConservesMomentum = {0:True,
                               1:False,
                               2:False,
                               3:True,
                               4:False}

#SHAKE error messages
ErrorMsgs = {10:"Constraint failure in RATTLE; large forces likely.",
             11:"Too many iterations in RATTLE (part 1)",
             12:"Too many iterations in RATTLE (part 2)",
             }

SourceVars = """
>>> inputmaps
int NSteps = self.StepsPerCycle

>>> outputmaps
int ErrorFlag = self.ErrorFlag

>>> defs
float dtsq2
float dt2
float idt
float Accel(Dim)
int i
int j
int istart
int istop
int ind
int Step
int ListInd
int m
bool Done
int MaxIter
int Iter
int MID
int AInd1
int AInd2
int a1
int a2
int b
float blsq
float iBoxL(Dim)
bool DoMinImage
bool Moving(MaxAtomsPerMol)
bool Moved(MaxAtomsPerMol)
float Pos0(NAtom,Dim)
float dPos1(Dim)
float dPos0(Dim)
float dVel(Dim)
float dfix(Dim)
float drSq1
float diffSq
float rpab
float gab
float rvab
float wc
float NH_wdti0
float NH_wdti1
float NH_wdti2 
float NH_wdti3
float NH_kT
float NH_NkT
float NH_scale
float AA
float NH_akin
int inos
float rn
float ranvec(Dim)
float dtfreq
float sqrtkt
float langevin1
float langevin2
float randomvelclip = -1.d0
bool ThisCalcVirial
bool ThisCalcDUParam
bool ThisCalcDWParam
int NActiveAxes
float Beta
float DeltaPos(Dim)
float CentroidPos(Dim)
float DeltaVol
float ScaleFactor(Dim)
float OldVol
float NewVol
float OldE
float NewE
float r
float lnP
float OldBoxL(Dim)
bool Acc
bool AxisMask(Dim)


>>> externals
ran2
ran2normarray
"""

Source_VVI = """
>>> Main
dtsq2 = TimeStep*TimeStep*0.5
dt2 = 0.5*TimeStep
idt = 1./TimeStep

do Step = 0, NSteps-1

    ThisCalcVirial = (CalcVirial .and. Step == NSteps - 1)
    ThisCalcDUParam = (CalcDUParam .and. Step == NSteps - 1)
    ThisCalcDWParam = (CalcDWParam .and. Step == NSteps - 1)
    
    [montecarlobarostat]

    [andersen]
    
    [removecom]
    
    [nosehoover]

    Pos0 = Pos                                                  || Pos0
    
    !velocity verlet integration part 1
    if (Thermostat == 4) then
        !langevin thermostat
        !based on Bussi and Parrinello, Physical Review E 75, 056707, 2007
        langevin1 = exp(-0.5d0 * LangevinGamma * TimeStep)
        langevin2 = sqrt((1.d0 - langevin1**2) * (kB * TempSet))
        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            Accel = Force(i,:) * iMass(i)
            call ran2normarray(Dim, ranvec)
            !clip the limits on the random variate for stability
            if (randomvelclip > 0.d0) then
                ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
            endif
            Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
            Pos(i,:) = Pos(i,:) + TimeStep*Vel(i,:) + dtsq2*Accel
            Vel(i,:) = Vel(i,:) + dt2*Accel
        enddo
    else
        !normal constant energy dynamics
        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            Accel = Force(i,:) * iMass(i)
            Pos(i,:) = Pos(i,:) + TimeStep*Vel(i,:) + dtsq2*Accel
            Vel(i,:) = Vel(i,:) + dt2*Accel
        enddo
    endif
    
    [rattle1]
    
    [external:calcenergyforces(Mode = 0, CalcForce = .true., CalcVirial = ThisCalcVirial, CalcDUParam = ThisCalcDUParam, CalcDWParam = ThisCalcDWParam)]
    
    !velocity verlet integration part 2
    
    if (Thermostat == 4) then
        !langevin thermostat
        langevin1 = exp(-0.5d0 * LangevinGamma * TimeStep)
        langevin2 = sqrt((1.d0 - langevin1**2) * (kB * TempSet))
        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            call ran2normarray(Dim, ranvec)
            !clip the limits on the random variate for stability
            if (randomvelclip > 0.d0) then
                ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
            endif
            Accel = Force(i,:) * iMass(i)
            Vel(i,:) = Vel(i,:) + dt2*Accel
            Vel(i,:) = langevin1*Vel(i,:) + langevin2*SqrtMass(i)*iMass(i)*ranvec
        enddo
    else
        !normal constant energy dynamics
        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            Accel = Force(i,:) * iMass(i)
            Vel(i,:) = Vel(i,:) + dt2*Accel
        enddo
    endif    
    
    [rattle2]
    
    [nosehoover]
    
    !update kinetic energy
    KEnergy = 0.d0
    do i = 0, NAtom-1
        if (.not. MolActive(MInd(i))==1) cycle
        KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
    enddo
    KEnergy = KEnergy * 0.5d0
    TEnergy = KEnergy + PEnergy
    
    !update running sums for total energy
    TEnergySum = TEnergySum + TEnergy
    TEnergySqSum = TEnergySqSum + TEnergy*TEnergy

enddo    


>>> montecarlobarostat
if (Barostat == 1 .and. mod(Step, BarostatStepFreq)==0) then
    Beta = 1. / (TempSet * kB)
    !update the attempt
    BarostatNAtt = BarostatNAtt + 1.
    
    NActiveAxes = count(BarostatUseAxis)
    OldBoxL = BoxL
    OldVol = product(BoxL)
    OldE = PEnergy
    [external:saveenergystate(Mode = 0)]
    
    !check if we need to scale independently
    if (.not. BarostatDoIsotropic) then
        !find a random active axis
        call ran2int(NActiveAxes, i)
        j = -1
        do ind = 0, Dim - 1
            if (BarostatUseAxis(ind)) j = j + 1
            if (j == i) exit
        enddo
        AxisMask = .false.
        AxisMask(ind) = .true.
    else
        AxisMask = BarostatUseAxis
    endif
    
    !choose a random volume change
    call ran2(r)
    DeltaVol = BarostatDeltaV * (2.d0 * r - 1.d0)
    NewVol = OldVol + DeltaVol
    
    if (NewVol > 0. .and. NewVol >= BarostatMinVol .and. NewVol <= BarostatMaxVol) then
       
        ScaleFactor = merge((NewVol / OldVol)**(1.d0 / dble(count(AxisMask))), 1.d0, AxisMask)
        BoxL = BoxL * ScaleFactor
        OldPos = Pos
            
        !now scale the molecule centers of mass
        do m = 0, NMol - 1
            !skip frozen and inactive mols
            if (MolActive(m) < 1) cycle
            
            !find the current centroid
            istart = MolRange(m)
            istop = MolRange(m+1) - 1
            CentroidPos = sum(Pos(istart:istop,:), dim=1) / AtomsPerMol(MolID(m))
            
            !find displacement
            DeltaPos = CentroidPos * (ScaleFactor - 1.d0)
            
            !update atom positions
            do i = istart, istop
                Pos(i,:) = Pos(i,:) + DeltaPos   
            enddo
        enddo
               
        !update energy 
        [external:calcenergyforces(Mode = 0, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)]
        NewE = PEnergy
        
        lnP = Beta * (OldE - NewE) - Beta * PresSet * DeltaVol
        lnP = lnP + NActiveMol * log(NewVol / OldVol)
        if (lnP >= 0) then
            Acc = .true.
        else
            call ran2(r)
            Acc = (exp(lnP) > r)
        endif
        
        if (Acc) then
            BarostatNAcc = BarostatNAcc + 1.
        else
            BoxL = OldBoxL
            Pos = OldPos
            [external:revertenergystate(Mode = 0)]
        endif     
    endif
endif

>>> andersen
!Andersen thermostats
if (Thermostat == 1) then
    !do an andersen massive collision update
    if (mod(AndersenStep, AndersenStepFreq)==0) then
        AndersenStep = 0
        sqrtkt = sqrt(kB * TempSet)
        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            call ran2normarray(Dim, ranvec)
            !clip the limits on the random variate for stability
            if (randomvelclip > 0.d0) then
                ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
            endif
            ranvec = ranvec * sqrtkt / sqrtMass(i)
            Vel(i,:) = ranvec
        enddo
    endif
    AndersenStep = AndersenStep + 1
endif

if (Thermostat == 2) then
    !do an andersen particle collision update
    dtfreq = TimeStep * AndersenCollisionFreq
    sqrtkt = sqrt(kB * TempSet)
    do m = 0, NMol-1
        if (.not. MolActive(m)==1) cycle
        call ran2(rn)
        if (rn < dtfreq) then
            do i = MolRange(m), MolRange(m+1)-1
                call ran2normarray(Dim, ranvec)
                !clip the limits on the random variate for stability
                if (randomvelclip > 0.d0) then
                    ranvec = merge(sign(randomvelclip, ranvec), ranvec, abs(ranvec) > randomvelclip)
                endif
                ranvec = ranvec * sqrtkt / sqrtMass(i)
                Vel(i,:) = ranvec
            enddo
        endif
    enddo
endif


>>> removecom
!remove the center of mass
if (RemoveCOMStepFreq > 0) then
    if (mod(RemoveCOMStep, RemoveCOMStepFreq)==0) then
        RemoveCOMStep = 0
        ranvec = 0.d0
        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            ranvec = ranvec + Mass(i) * Vel(i,:)
        enddo
        ranvec = ranvec / NAtom
        do i = 0, NAtom-1
            if (.not. MolActive(MInd(i))==1) cycle
            Vel(i,:) = Vel(i,:) - ranvec * iMass(i)
        enddo
    endif
    RemoveCOMStep = RemoveCOMStep + 1
endif

    
>>> nosehoover
!NOSE-HOOVER ROUTINES
if (Thermostat == 3) then
    !update kinetic energy
    KEnergy = 0.d0
    do i = 0, NAtom-1
        KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
    enddo
    KEnergy = KEnergy * 0.5d0
    TEnergy = KEnergy + PEnergy
    !set frequently used variables
    NH_wdti0 = TimeStep
    NH_wdti1 = TimeStep * 0.5d0
    NH_wdti2 = TimeStep * 0.25d0
    NH_wdti3 = TimeStep * 0.125
    NH_kT = TempSet * kB 
    NH_NkT = NH_kT * dble(NDOF - Dim)
    NH_scale = 1.D0  
    !get kinetic energy
    NH_akin = 2.d0 * KEnergy
    !update the forces
    Glogs(0) = (NH_akin - NH_NkT) / QMass(0)
    !update the thermostat velocities
    Vlogs(NH_N-1) = Vlogs(NH_N-1) + Glogs(NH_N-1) * NH_wdti2
    do inos = 1, NH_N - 1
        AA = exp( -NH_wdti3 * Vlogs(NH_N-inos) )
        Vlogs(NH_N-inos-1) = Vlogs(NH_N-inos-1)*AA*AA + NH_wdti2*Glogs(NH_N-inos-1)*AA
    enddo
    !update the particle velocities
    AA = exp( -NH_wdti1 * Vlogs(0) )
    NH_scale = NH_scale * AA
    !update the forces
    Glogs(0) = (NH_scale*NH_scale*NH_akin - NH_NkT) / QMass(0)
    !update the thermostat positions
    do inos = 0, NH_N - 1
        Xlogs(inos) = Xlogs(inos) + Vlogs(inos) * NH_wdti1
    enddo
    !update the thermostat velocities
    do inos = 1, NH_N-1
        AA = exp( -NH_wdti3 * Vlogs(inos) )
        Vlogs(inos-1) = Vlogs(inos-1)*AA*AA + NH_wdti2*Glogs(inos-1)*AA
        Glogs(inos) = (QMass(inos-1)*Vlogs(inos-1)*Vlogs(inos-1)-NH_kT) / QMass(inos)
    enddo
    Vlogs(NH_N-1) = Vlogs(NH_N-1) + NH_wdti2*Glogs(NH_N-1)
    !update the particle velocities
    if (NH_scale > 0.) then
        Vel = Vel * NH_scale
    endif
endif    
"""

Source_RattleI = """
>>> rattle1
!rattle 1
DoMinImage = any(BoxL > 0.d0)
Done = .false.
iBoxL = 1.d0 / max(1.d-300, BoxL)

!loop over molecules
do m = 0, NMol-1

    MID = MolID(m)
    
    !check if there are no rigid bonds
    if (RBondRange(MID+1) - RBondRange(MID) == 0) cycle
    
    !compute max iterations
    MaxIter = RattleMaxIterPerAtom * (MolRange(m+1) - MolRange(m))  
    
    Iter = 0
    Done = .false.
    Moving = .false.
    Moved = .true.

    do while ((.not. Done) .and. Iter < MaxIter) 
        Done = .true.

        do b = RBondRange(MID), RBondRange(MID+1) - 1
            AInd1 = RBondInd(b,0)
            AInd2 = RBondInd(b,1)
            a1 = AInd1 + MolRange(m)
            a2 = AInd2 + MolRange(m)
            blsq = RBondLengthSq(b)

            if (Moved(AInd1) .or. Moved(AInd2) ) then
                dPos1 = Pos(a1,:) - Pos(a2,:)
                if (DoMinImage) dPos1 = dPos1 - dnint(dPos1 * iBoxL) * BoxL
                drSq1 = sum(dPos1*dPos1)
                diffsq = blsq - drSq1

                if (abs(diffsq) > (blsq * RattleTol1)) then
                    dPos0 = Pos0(a1,:) - Pos0(a2,:)
                    if (DoMinImage) dPos0 = dPos0 - dnint(dPos0 * iBoxL) * BoxL
                    rpab = dot_product(dPos1,dPos0)
                    if (rpab < blsq * RattleTol2) then
                        !constraint failure
                        ErrorFlag = 10
                        return
                    endif
                    gab = diffsq / ( 2.d0 * (iMass(a1) + iMass(a2)) * rpab )
                    dfix = dPos0 * gab
                    Pos(a1,:) = Pos(a1,:) + iMass(a1) * dfix
                    Pos(a2,:) = Pos(a2,:) - iMass(a2) * dfix
                    dfix = dfix * idt
                    Vel(a1,:) = Vel(a1,:) + iMass(a1) * dfix
                    Vel(a2,:) = Vel(a2,:) - iMass(a2) * dfix
                    Moving(AInd1) = .true.
                    Moving(AInd2) = .true.
                    Done = .false.
                endif
            endif
        enddo
        Moved = Moving
        Moving = .false.
        Iter = Iter + 1
    enddo
    
    !end of iterative loop
    if (.not. Done) then
        !too many iterations
        ErrorFlag = 11
        return
    endif

enddo


>>> rattle2
!rattle 2
Done = .false.
wc = 0.d0
iBoxL = 1.d0 / max(1.d-300, BoxL)

!loop over molecules and bonds
do m = 0, NMol - 1

    MID = MolID(m)
    
    !check if there are no rigid bonds
    if (RBondRange(MID+1) - RBondRange(MID) == 0) cycle
    
    !compute max iterations
    MaxIter = RattleMaxIterPerAtom * (MolRange(m+1) - MolRange(m))  
    
    Iter = 0
    Done = .false.
    Moving = .false.
    Moved = .true.

    do while ((.not. Done) .and. Iter < MaxIter) 
        Done = .true.

        do b = RBondRange(MID), RBondRange(MID+1) - 1
            AInd1 = RBondInd(b,0)
            AInd2 = RBondInd(b,1)
            a1 = AInd1 + MolRange(m)
            a2 = AInd2 + MolRange(m)
            blsq = RBondLengthSq(b)

            if (Moved(AInd1) .or. Moved(AInd2) ) then
                dVel = Vel(a1,:) - Vel(a2,:)
                dPos1 = Pos(a1,:) - Pos(a2,:)
                if (DoMinImage) dPos1 = dPos1 - anint(dPos1 * iBoxL) * BoxL
                rvab = dot_product(dPos1, dVel)

                if (abs(rvab*TimeStep) > (RattleTol3 * blsq)) then
                    gab = -rvab / ((iMass(a1) + iMass(a2)) * blsq)
                    wc = wc + gab * blsq
                    dfix = dPos1 * gab
                    Vel(a1,:) = Vel(a1,:) + iMass(a1) * dfix
                    Vel(a2,:) = Vel(a2,:) - iMass(a2) * dfix
                    Moving(AInd1) = .true.
                    Moving(AInd2) = .true.
                    Done = .false.
                endif
            endif
        enddo
        Moved = Moving
        Moving = .false.
        Iter = Iter + 1
    enddo
    
    !end of iterative loop
    if (.not. Done) then
        !too many iterations
        ErrorFlag = 12
        return
    endif
enddo

!update the virial
wc = wc / dt2
Virial = Virial - wc
"""


Source_NoRattle = """
>>> rattle1
!no rattle for this system

>>> rattle2
!no rattle for this system
"""


Source_VVQ = """
>>> Main
!velocity verlet quench 1
dtsq2 = TimeStep*TimeStep*0.5
dt2 = 0.5*TimeStep
idt = 1./TimeStep

do Step = 0, NSteps-1
    Pos0 = Pos                                                  || Pos0
    KEnergy = 0.d0
    TEnergy = KEnergy + PEnergy
    ThisCalcVirial= (CalcVirial .and. Step == NSteps - 1)
    ThisCalcDUParam = (CalcDUParam .and. Step == NSteps - 1)
    ThisCalcDWParam = (CalcDWParam .and. Step == NSteps - 1)
    
    do i = 0, NAtom-1
        if (.not. MolActive(MInd(i))==1) cycle
        Accel = Force(i,:) * iMass(i)
        Vel(i,:) = 0.d0
        Pos(i,:) = Pos(i,:) + dtsq2*Accel
    enddo
       
    [rattle1]
    
    [external:calcenergyforces(Mode = 0, CalcForce = .true., CalcVirial = ThisCalcVirial, CalcDUParam = ThisCalcDUParam, CalcDWParam = ThisCalcDWParam)]
enddo
"""

Source_RattleQ = """
>>> rattle1
!rattle quench 1
Done = .false.
iBoxL = 1.d0 / max(1.d-300, BoxL)

!loop over molecules
do m = 0, NMol-1

    MID = MolID(m)
    
    !check if there are no rigid bonds
    if (RBondRange(MID+1) - RBondRange(MID) == 0) cycle
    
    !compute max iterations
    MaxIter = RattleMaxIterPerAtom * (MolRange(m+1) - MolRange(m))  
    
    Iter = 0
    Done = .false.
    Moving = .false.
    Moved = .true.

    do while ((.not. Done) .and. Iter < MaxIter) 
        Done = .true.

        do b = RBondRange(MID), RBondRange(MID+1) - 1
            AInd1 = RBondInd(b,0)
            AInd2 = RBondInd(b,1)
            a1 = AInd1 + MolRange(m)
            a2 = AInd2 + MolRange(m)
            blsq = RBondLengthSq(b)

            if (Moved(AInd1) .or. Moved(AInd2) ) then
                dPos1 = Pos(a1,:) - Pos(a2,:)
                if (DoMinImage) dPos1 = dPos1 - dnint (dPos1 * iBoxL) * BoxL
                drSq1 = sum(dPos1*dPos1)
                diffsq = blsq - drSq1
                
                if (abs(diffsq) > (blsq * RattleTol1)) then
                    dPos0 = Pos0(a1,:) - Pos0(a2,:)
                    if (DoMinImage) dPos0 = dPos0 - dnint (dPos0 * iBoxL) * BoxL
                    rpab = dot_product(dPos1,dPos0)
                    if (rpab < blsq * RattleTol2) then
                        !constraint failure
                        ErrorFlag = 10
                        return
                    endif
                    gab = diffsq / ( 2.d0 * (iMass(a1) + iMass(a2)) * rpab )
                    dfix = dPos0 * gab
                    Pos(a1,:) = Pos(a1,:) + iMass(a1) * dfix
                    Pos(a2,:) = Pos(a2,:) - iMass(a2) * dfix
                    Moving(AInd1) = .true.
                    Moving(AInd2) = .true.
                    Done = .false.
                endif
            endif
        enddo
        Moved = Moving
        Moving = .false.
        Iter = Iter + 1
    enddo
    
    !end of iterative loop
    if (.not. Done) then
        !too many iterations
        ErrorFlag = 11
        return
    endif

enddo
"""

Source_KE = """
!update kinetic energy
KEnergy = 0.d0
do i = 0, NAtom-1
    if (.not. MolActive(MInd(i))==1) cycle
    KEnergy = KEnergy + dot_product(Vel(i,:), Vel(i,:)) * Mass(i)
enddo
KEnergy = KEnergy * 0.5d0
TEnergy = KEnergy + PEnergy
"""


def GetKEFortCode(Sys, LibVars0):
    from sim.fortran import ReplaceToken
    Source = SourceVars + Source_KE
    for (Var, ModuleVar) in LibVars0.items():
        Source = ReplaceToken(Source, Var, ModuleVar)
    fc = FortCode(Source)
    return fc  

def GetVVIFortCode(Sys, LibVars0):
    from sim.fortran import ReplaceToken
    HasRigid = any([MolType.HasRigid for MolType in Sys.World])
    if HasRigid:
        Source = SourceVars + Source_VVI + Source_RattleI
    else:
        Source = SourceVars + Source_VVI + Source_NoRattle
    for (Var, ModuleVar) in LibVars0.items():
        Source = ReplaceToken(Source, Var, ModuleVar)
    fc = FortCode(Source)
    return fc  

def GetVVQFortCode(Sys, LibVars0):
    from sim.fortran import ReplaceToken
    HasRigid = any([MolType.HasRigid for MolType in Sys.World])
    if HasRigid:
        Source = SourceVars + Source_VVQ + Source_RattleQ
    else:
        Source = SourceVars + Source_VVQ + Source_NoRattle
    for (Var, ModuleVar) in LibVars0.items():
        Source = ReplaceToken(Source, Var, ModuleVar)
    fc = FortCode(Source)
    return fc  
          




