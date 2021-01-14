#/usr/bin/env python


### Translates Sys objects to LAMMPS input files, in SIM suite.
### coded by SPC, MSS

import os

import numpy as np

import sim.potential.base.potentialtypes as ptypes

OMP_NumThread = 10
#executable
LammpsExec = os.environ.get("LAMMPSEXEC", "lmp")
UseOMP = True
OMP_NumThread=2

#debug options
WARNINGS = True


#defaults for this class
NBondPotentialBins = 500
NPairPotentialBins = 500
NAnglePotentialBins = 500
NDihedralPotentialBins = 500
NLocalDensityPotentialBins = 500
NOtherPotentialBins = 500
TableInterpolationStyle = 'spline'

#enable for spline bonded interaction
BondSpline = False

#pair spline inner cutoff in Angstroms
InnerCutoff = 0.1

#use table2 pair style which allows for soft potentials?
UseTable2 = False

#barostat damping factor
PDampFactor = 1000.
#neighbor list settings
NeighDelay = None
NeighEvery = None
NeighCheck = None
SkinDistance = None
NeighPage = None
NeighOne = None #2000 is default
#max pair energy scale in units of kB T
MaxPairEnekBT = 20.

#conversion factor for real (kcal/mol) to metal (eV) units
MetalConv = 0.04336210 

#reset positions before running Lammps
ResetPosFirst = False

#Long range electrostatics option
KspaceStyleCmd = "pppm 1.e-4"

#Shake parameters
ShakeTol = 1.e-4
ShakeIter = 50
#rigid bond/angle force constant (used during minimization); units are EScale/Angstrom^2
RBondFConst = 100.
RAngFConst = 100.





#make an error class 
class LammpsError(Exception):
  def __init__(self, Msg):
    self.Msg = Msg
  def __str__(self):
    return str(self.Msg)
    
 

def AName(Atom1):
    """Returns a string name of an atom."""
    #get the atom name and append the index in the molecule
    return Atom1.Name.strip() + "%d" % (Atom1.AInd + 1)
    
    
def PName(Sys, P):
    """Returns a string name of a potential."""
    if type(P) is list:
        PList = P
        s = "PairPotential"
        for (P, TypeInd) in PList:
            s += "%d[%d]" % (Sys.ForceField.index(P), TypeInd)
        return s
    else:
        i = Sys.ForceField.index(P) + 1
        if P.Type == ptypes.PairPotential:
            if P.Filter.Bonded:
                return "BondPotential%d" % i
            else:
                return "PairPotential%d" % i
        elif P.Type == ptypes.AnglePotential:
            return "AnglePotential%d" % i
        elif P.Type == ptypes.TorsionPotential:
            return "DihedralPotential%d" % i
        elif P.Type == ptypes.LocalDensityPotential:
            return "LocalRhoPotential%d" % i
        else:
            raise LammpsError("Don't recognize potential type.")
            
def PotentialTableBond(Sys, P, Cut):
    """Returns a string of potential arguments and values."""
    #start the string; #make a label
    s = "\n%s\n" % PName(Sys, P)
    N = NBondPotentialBins
    #bond potential
    #HERE RAISE FLAG IF NO CUTOFF
    if Cut is None:
        xmax = np.max([P.Arg.ReportMax.max() for (P, TypeInd) in PList])
    else:
        xmax = Cut
    xmin = 0.
    #if UseTable2:
    #    xmin = 0.
    #else:
    #    xmin = Sys.Units.Angstrom * InnerCutoff
    xvals = np.linspace(xmin, xmax, N)
    uvals = np.zeros_like(xvals)
    duvals = np.zeros_like(xvals)

    uvals += np.array([P.Val(x) for x in xvals])
    duvals += np.array([-P.DVal(x) for x in xvals])

    if not np.isfinite(uvals[0]): uvals[0] = 1.e300

    #get the maximum pair energy
    #MaxPairEne = Sys.Units.kB * Sys.TempSet * MaxPairEnekBT
    #indices where energy is greater
    #ind = np.where(uvals > MaxPairEne)[0]
    #if len(ind):
        #find the first index where energy is valid
    #    i = ind[-1] + 1
        #do a linear extrapolation in the hard core region
    #    uvals[:i] = (xvals[i] - xvals[:i]) * duvals[i] + uvals[i]
    #    duvals[:i] = duvals[i]
    s += "N %d \n\n" % (N)
    for (j, (x, u, du)) in enumerate(zip(xvals, uvals, duvals)):
        s += "%-4d %-12.5e %-12.5e %12.5e\n" % (j+1, x, u, du)
    return s
        

def PotentialTablePair(Sys, PList, Cut):
    """Returns a string of potential arguments and values."""  
    #start the string; #make a label
    s = "\n%s\n" % PName(Sys, PList)  
    N = NPairPotentialBins
    #pair potential
    #HERE RAISE FLAG IF NO CUTOFF
    if Cut is None:
        xmax = np.max([P.Arg.ReportMax.max() for (P, TypeInd) in PList])
    else:
        xmax = Cut
    if UseTable2:
        xmin = 0.
    else:
        xmin = Sys.Units.Angstrom * InnerCutoff
    xvals = np.linspace(xmin, xmax, N)
    uvals = np.zeros_like(xvals)
    duvals = np.zeros_like(xvals)
    for (P, TypeInd) in PList:
        P.SetTypeInd(TypeInd)
        uvals += np.array([P.Val(x) for x in xvals])
        duvals += np.array([-P.DVal(x) for x in xvals])
    if not np.isfinite(uvals[0]): uvals[0] = 1.e300
    #get the maximum pair energy
    MaxPairEne = Sys.Units.kB * Sys.TempSet * MaxPairEnekBT
    #indices where energy is greater
    ind = np.where(uvals > MaxPairEne)[0]
    if len(ind):
        #find the first index where energy is valid
        i = ind[-1] + 1
        #do a linear extrapolation in the hard core region
        uvals[:i] = (xvals[i] - xvals[:i]) * duvals[i] + uvals[i]
        duvals[:i] = duvals[i]
    s += "N %d R %12.5e %12.5e\n\n" % (N, xmin, xmax)
    for (j, (x, u, du)) in enumerate(zip(xvals, uvals, duvals)):
        s += "%-4d %-12.5e %-12.5e %12.5e\n" % (j+1, x, u, du) 
    return s


def PotentialTableNullPair(Cut, fmt = r"%-11.4e"):
    if Cut is None: Cut = 10.
    N = NPairPotentialBins
    s = '\nNonBondNull\nN %d R 0.0000000001 %12.5e\n\n' % (N, Cut)
    bins = np.linspace(0., 10., N)
    for (j, x) in enumerate(bins):
        s += "%-4d" % (j+1) + fmt % x + " " + fmt % 0.0 + " " + fmt % 0.0 + "\n"
    return s
            

def PotentialTableAngle(Sys, P):
    """Returns a string of potential arguments and values."""     
    #start the string; #make a label
    s = "\n%s\n" % PName(Sys, P)      
    #make the table 
    N = NAnglePotentialBins
    #angle potential
    xscale = 180. / np.pi
    xvals = np.linspace(0, np.pi, N)
    uvals = np.array([P.Val(x) for x in xvals])
    duvals = np.array([-P.DVal(x) / xscale for x in xvals])
    xvals = xvals * xscale
    s += "N %d \n\n" % N
    for (j, (x, u, du)) in enumerate(zip(xvals, uvals, duvals)):
        s += "%-4d %-12.5e %-12.5e %12.5e\n" % (j+1, x, u, du) 
    return s

def PotentialTableAngleRigid(Sys, Name, Theta0):
    """Returns a string of potential arguments and values."""     
    #start the string; #make a label
    s = "\n%s\n" % Name      
    #make the table 
    N = NAnglePotentialBins
    #force constant
    k = RAngFConst * Sys.Units.EScale
    #angle potential
    xvals = np.linspace(0., 180., N)
    #do a harmonic potential
    uvals = np.array([0.5*k*(x-Theta0)**2 for x in xvals])
    duvals = np.array([-0.5*k*(x-Theta0) for x in xvals])
    s += "N %d EQ %12.5e\n\n" % (N, Theta0 * np.pi / 180.)
    for (j, (x, u, du)) in enumerate(zip(xvals, uvals, duvals)):
        s += "%-4d %-12.5e %-12.5e %12.5e\n" % (j+1, x, u, du) 
    return s
    
    
def PotentialTableDihedral(Sys, P):
    """Returns a string of potential arguments and values.""" 
    #start the string; #make a label
    s = "\n\n%s\n" % PName(Sys, P)      
    #make the table 
    N = NDihedralPotentialBins        
    #dihedral potential
    #The explicit pi below ensures the maximum angle is rounded down from 2*pi.
    #Otherwise LAMMPS complains that the maximum angle is too big when string formatting rounds up.
    xvals = np.linspace(0., 6.28318, N)
    uvals = np.array([P.Val(x) for x in xvals])
    s += "N %d RADIANS NOF\n" % N
    for (j, (x, u)) in enumerate(zip(xvals, uvals)):
        s += "%-4d %-12.5e %-12.5e\n" % (j+1, x, u) 
    return s 

    

def MakeLocalDensityFile(LocalDensityFile, LD):
    """Uses LAMMPS custom localdensity potential to do local density potentials."""
    
    N = NLocalDensityPotentialBins
        
    s = "#LOCAL DENSITY POTENTIALS\n\n"
    
    #add number of ld potentials and number of bins
    s += "%d %d\n" % (len(LD), N)
    s += "\n"
    
    #output local density potentials
    for (ind, (P, CenterIDs, NeighIDs)) in enumerate(LD):
        #add lower and upper cutoffs
        s += "%14.7e %14.7e\n" % (P.InnerCut[0], P.Cut)
        #add central atom types
        s += " ".join([str(ID+1) for ID in CenterIDs]) + "\n"
        #add neighbor atom types
        s += " ".join([str(ID+1) for ID in NeighIDs]) + "\n"
        #add rhomin rhomax deltarho
        dRho = (P.RhoMax - P.RhoMin) / (N - 1)
        s += "%14.7e %14.7e %14.7e\n" % (P.RhoMin, P.RhoMax, dRho)
        #add the potential values F(rho)
        Frho = [P.Val(P.RhoMin + i*dRho) for i in xrange(N)]
        for Frhoi in Frho:
            s += "%14.7e\n" % Frhoi
        s += "\n"
    
    #write out the file
    file(LocalDensityFile, "w").write(s)
    


def CheckForLammps(Sys):
    """Checks a system to make sure that it can be prepped for LAMMPS.
The Has varaiables are True/False for potential types.
The others are special bonds coefficients for 1-2, 1-3, and 1-4 interactions."""
    import sim.potential
    import sim.atomselect    
        
    HasBond = False
    HasAngle = False
    HasDihedral = False
    HasNonbond = False
    HasLocalRho = False
    HasEwald = False
    HasRigid = False
    HasExtField = False
    
    #check force field types
    UseSites = False
    for P in Sys.ForceField:
        if P.Filter.HasSiteType():
            UseSites = True
        if P.Arg.NType > 1 and (P.Type != ptypes.PairPotential or P.Filter.Bonded):
            raise LammpsError("Cannot convert multi-type potentials besides pair nonbonded.")
        if isinstance(P, sim.potential.Ewald) or isinstance(P, sim.potential.ewald.Ewald):
            HasEwald = True
        elif P.Type == ptypes.PairPotential and P.Filter.Bonded:
            if isinstance(P, sim.potential.Bond) or isinstance(P, sim.potential.PairSpline):
                HasBond = True
            else:
                raise LammpsError("Cannot convert bonded potentials that are not harmonic or spline.")
        elif P.Type == ptypes.PairPotential:
            HasNonbond = True
        elif P.Type == ptypes.AnglePotential:
            HasAngle = True
        elif P.Type == ptypes.TorsionPotential:
            HasDihedral = True
        elif P.Type == ptypes.LocalDensityPotential:
            HasLocalRho = True
        elif P.Type == ptypes.FieldPotential:
            if P.Names[0] == "external_sinusoid":
                HasExtField = True
            else:
                raise LammpsError("External field only implemented for sinusoid so far")
        else:
            #NOTE: if we add field potentials to this, then need to change recentering of positions
            raise LammpsError("Can't convert potential %s" % P.Name)
    print('HasEwald {}'.format(HasEwald))
    #check for rigid or constraints
    for Mol in Sys.Mol:
        if Mol.HasRigid:
            HasRigid = True
            HasBond = True 
            HasAngle = True
            break
            
    #calculate the number of types
    if UseSites:
        Types = Sys.World.SiteTypes
        def GetID(Type1):
            return Type1.SID
        if WARNINGS: 
            print "WARNING: Using site rather than atom types for LAMMPS."
            print "         LAMMPS performance may be slow."
    else: 
        Types = Sys.World.AtomTypes
        def GetID(Type1):
            return Type1.AID
            
    #make sure that all bonded site types have a harmonic potential or a rigid bond
    #first build a bond matrix from potentials
    STypes = Sys.World.SiteTypes
    NSTypes = len(STypes)
    BondMatrix = np.zeros((NSTypes, NSTypes), dtype = bool)
    #first mark the rigid bonds
    BondOrdMatrix = Sys.World.GetBondOrdMatrix(ShowRigid=True)
    BondMatrix[BondOrdMatrix == -2] = True
    #now find the ones with a potential
    for P in Sys.ForceField:
        if P.Type == ptypes.PairPotential and P.Filter.Bonded:
            for (SType1, SType2) in P.Filter.Select(STypes):
                BondMatrix[SType1.SID, SType2.SID] = True
                BondMatrix[SType2.SID, SType1.SID] = True
    #now run through all bonds and ensure each had a potential
    for (SType1, SType2) in sim.atomselect.BondPairs.Select(STypes):
        ID1, ID2 = SType1.SID, SType2.SID
        if not BondMatrix[ID1, ID2]:
            raise LammpsError("Could not find a harmonic or rigid bond between site types %d and %d" % (ID1, ID2))

    def AppliesToIntramolecular(Filter):
        for Mol in Sys.Mol:
            Atoms = [a for a in Mol]
            for Group in Filter.Select(Atoms):
                #made it into the loop; there is at least one group
                return True
        return False
            
    #find a global cutoff -- this must be used if kspace style is used
    GlobalCut = None
    #check cutoff lengths
    if np.any(Sys.BoxL > 0):
        MaxCut = 0.5 * np.min(Sys.BoxL[Sys.BoxL > 0])
    else:
        MaxCut = 1.e300
    if HasEwald:
        for P in Sys.ForceField:
            #check local density
            if not P.Type == ptypes.PairPotential or P.Filter.Bonded: continue
            #check that the cutoff isn't too big
            if not P.Cut is None and P.Cut > MaxCut:
                raise LammpsError("Cutoff for potential %s greater (%11.4e) than half min box length (%11.4e)." % (P.Name, P.Cut, MaxCut))
            #check the cutoff against global
            if GlobalCut is None:
                GlobalCut = P.Cut
            elif not GlobalCut == P.Cut:
                raise LammpsError("When using Ewald, all pair potentials must have same cutoff as real space.")
    
    #determine how to apply nonbonded interactions to 1-2, 1-3, and 1-4 pairs
    LJCoef12, LJCoef13, LJCoef14 = None, None, None
    for P in Sys.ForceField:
        #skip Ewald -- has separate special bonds in LAMMPS
        if isinstance(P, sim.potential.Ewald): continue
    
        #check filters on atoms in system
        F = P.Filter

        #check local density
        if P.Type == ptypes.LocalDensityPotential:
            if F.HasIntra():
                raise LammpsError("Cannot process potentials with Intra = True or False (%s)." % P.Name)
                
        #skip if not a pair potential
        if not P.Type == ptypes.PairPotential or P.Filter.Bonded: continue

        #skip if there are no relevant intramolecular atom pairs
        if not AppliesToIntramolecular(F): continue
    
        #check if there is an intra or intermolecular specification; LAMMPS cannot do
        if F.HasIntra():
            raise LammpsError("Cannot process potentials with Intra = True or False (%s)." % P.Name)
               
        #get the bond order range for this potential
        if F.MinBondOrd is None:
            MinBondOrd = 0
        else:
            if F.MinBondOrd > 5:
                raise LammpsError("Cannot process potentials with minimum bond order > 4 (%s)" % P.Name)
            MinBondOrd = F.MinBondOrd
        if not F.MaxBondOrd is None:
            raise LammpsError("Cannot process potentials with a maximum bond order set (%s)" % P.Name)
        #check usage of special interactions
        ThisLJCoef12 = float(2 >= MinBondOrd)
        ThisLJCoef13 = float(3 >= MinBondOrd)
        ThisLJCoef14 = float(4 >= MinBondOrd)
        #now see how to set coefficients; check that we haven't already set since
        #these are global in lammps
        if LJCoef12 is None:
            LJCoef12 = ThisLJCoef12
        elif not LJCoef12 == ThisLJCoef12:
            raise ValueError("Found two nonbonded potentials with different special 1-2 bond treatments.")
        if LJCoef13 is None:
            LJCoef13 = ThisLJCoef13
        elif not LJCoef13 == ThisLJCoef13:
            raise ValueError("Found two nonbonded potentials with different special 1-3 bond treatments.")
        if LJCoef14 is None:
            LJCoef14 = ThisLJCoef14
        elif not LJCoef14 == ThisLJCoef14:
            raise ValueError("Found two nonbonded potentials with different special 1-4 bond treatments.")
    if LJCoef12 is None: LJCoef12 = 1.
    if LJCoef13 is None: LJCoef13 = 1.
    if LJCoef14 is None: LJCoef14 = 1.
            
    #now check ewald coefs
    EWCoef12, EWCoef13, EWCoef14 = 1., 1., 1.
    for P in Sys.ForceField:
        if isinstance(P, sim.potential.Ewald):
            if P.ExcludeBondOrd >= 2:
                EWCoef12 = 0.
            if P.ExcludeBondOrd >= 3:
                EWCoef13 = 0.
            if P.ExcludeBondOrd >= 4:
                EWCoef14 = 0.
            if P.ExcludeBondOrd > 4:
                raise LammpsError("Cannot process potentials with minimum bond order > 4 (%s)" % P.Name)

    #return 
    return (HasBond, HasAngle, HasDihedral, HasNonbond, HasLocalRho, HasEwald, HasRigid, HasExtField,
            GlobalCut, LJCoef12, LJCoef13, LJCoef14, EWCoef12, EWCoef13, EWCoef14, Types, GetID)


def AtomsIndexKey(Atoms):
    key = [a.Ind for a in Atoms]
    if len(key) > 1:
        if key[-1] < key[1]:
            key = reversed(key)
    return tuple(key)
    

def MakeLammps(Sys, Prefix = "", InFile = "lammps.in", DataFile = "lammps.data", RestartFile = None,
               PairFile = "lammps.pair", AngleFile = "lammps.angle", DihedralFile = "lammps.dihedral",
               LocalDensityFile = "lammps.localdensity", BondFile = "lammps.bond", LammpsCommands = "", LammpsHeader = ""):   
    """Makes LAMMPS input files for Sys.
Returns the filenames used: InFile, DataFile, PairFile, AngleFile, DihedralFile"""
    import sim.units
    
    #check if we need to add file prefixes to this run
    InFile = Prefix + InFile
    DataFile = Prefix + DataFile
    PairFile = Prefix + PairFile
    AngleFile = Prefix + AngleFile  
    DihedralFile = Prefix + DihedralFile  
    BondFile = Prefix + BondFile 
    LocalDensityFile = Prefix + LocalDensityFile      
       
    #check and get information for preppring lammps        
    (HasBond, HasAngle, HasDihedral, HasNonbond, HasLocalRho, HasEwald, HasRigid, HasExtField, 
     GlobalCut, LJCoef12, LJCoef13, LJCoef14, EWCoef12, EWCoef13, EWCoef14, 
     Types, GetID) = CheckForLammps(Sys) 
    NTypes = len(Types)
    STypes = Sys.World.SiteTypes

    #Build LAMMPS input file
    #It turns out this will need to be a stand alone molecule data file 
    #The file will contain all of the information about bond topology etc.      
    
    #choose a random number seed    
    RanSeed = 100000 * np.random.rand()

    #choose boundary type
    Boundary = ''
    for L in Sys.BoxL:
        if L > 0:
            Boundary += 'p '
        else:
            Boundary += 's '

    #choose the units type and time scale
    if Sys.Units == sim.units.DimensionlessUnits:
        Units = 'lj'
    elif Sys.Units == sim.units.MKSUnits:
        Units = 'si'
    elif Sys.Units == sim.units.AtomicUnits:
        Units = 'real'
    else:
        raise LammpsError("Don't recognize type of units in system.")
        
    #make atom names
    AtomNames = " ".join([Type1.Name.strip() for Type1 in Types])
    
    #atom style
    if HasEwald:
        AtomStyle = "full"
    else:
        AtomStyle = "molecular"
    
    #dictionary of tokens
    d = {"BOUNDARY" : Boundary, 
         "UNITS" : Units,
         "DATAFILE" : DataFile, 
         "RANSEED" : RanSeed,
         "NBONDPOTENTIALBINS" : NBondPotentialBins,
         "NPAIRPOTENTIALBINS" : NPairPotentialBins,
         "NANGLEPOTENTIALBINS" : NAnglePotentialBins,
         "NDIHEDRALPOTENTIALBINS" : NDihedralPotentialBins,
         "TEMP" : Sys.TempSet,
         "ATOMNAMES" : AtomNames,
         "ATOMSTYLE" : AtomStyle,
         "GLOBALCUT" : GlobalCut,
         "KSPACESTYLE" : KspaceStyleCmd,
         "TABLEINTERPOLATION" : TableInterpolationStyle,
         "RESTARTFILE": RestartFile}
    
    #add the table style
    if UseTable2:
        d["TABLESTYLE"] = "table2"
    else:
        d["TABLESTYLE"] = "table"

    #start input string
    if len(LammpsHeader):
        s = LammpsHeader + "\n"
    else:
        s = ""
        
    s += """# LAMMPS INPUT FILE GENERATED BY SIM SUITE
clear

# Initialize simulation box
"""
    if UseOMP:
        s += """package omp     {}\n""".format(OMP_NumThread)
        
    s += """dimension       3
boundary        %(BOUNDARY)s
units           %(UNITS)s
atom_style      %(ATOMSTYLE)s

# Set potential styles
""" % d
    if HasLocalRho:
        if HasEwald:
            s += "pair_style      hybrid/overlay %(TABLESTYLE)s spline %(NPAIRPOTENTIALBINS)d coul/long %(GLOBALCUT)12.5e localdensity\n" % d
        else:
            s += "pair_style      hybrid/overlay %(TABLESTYLE)s spline %(NPAIRPOTENTIALBINS)d localdensity\n" % d
    else:
        if HasEwald:
            s += "pair_style      hybrid/overlay %(TABLESTYLE)s spline %(NPAIRPOTENTIALBINS)d coul/long %(GLOBALCUT)12.5e\n" % d
        else:
            s += "pair_style      %(TABLESTYLE)s %(TABLEINTERPOLATION)s %(NPAIRPOTENTIALBINS)d\n" % d
    if HasEwald:
        s += "kspace_style    %(KSPACESTYLE)s\n" % d
    if HasBond:
        if BondSpline:
            s += "bond_style      table spline %(NBONDPOTENTIALBINS)d\n" % d
        else: 
            s += "bond_style      harmonic\n" % d
    if HasAngle:
        s += "angle_style     table spline %(NANGLEPOTENTIALBINS)d\n" % d
    if HasDihedral:
        s += "dihedral_style  table spline %(NDIHEDRALPOTENTIALBINS)d\n" % d
    
    if not RestartFile == None:
        s +="""
# Read molecule data and set initial velocities
read_restart       %(RESTARTFILE)s 

# Assign potentials
""" % d
    else:
        s +="""
# Read molecule data and set initial velocities
read_data       %(DATAFILE)s 
velocity        all create %(TEMP)11.4e %(RANSEED)d rot yes dist gaussian

# Assign potentials
""" % d

    #add ewald pair coeff
    if HasEwald:
        s += "pair_coeff      *     *     coul/long\n"
        
    #find all nonbonded interactions
    #make a matrix of indices to potentials
    PIndex = [] #2D array indexed by atom types; for each atom pair give list of pair potentials (and an index pointing to the pair) involving said atom pair
    for ID1 in range(NTypes):
        PIndex.append([])
        for ID2 in range(0, ID1+1):
            PIndex[-1].append([])
    #loop through nonbond pair potentials
    for (i, P) in enumerate(Sys.ForceField):
        if isinstance(P, sim.potential.Ewald): continue
        if P.Type == ptypes.PairPotential and not P.Filter.Bonded:
            #select atom types that work for this potential
            for (SType1, SType2) in P.Filter.Select(STypes):
                ID1, ID2 = GetID(SType1), GetID(SType2)
                if ID1 < ID2:
                    ID1, ID2 = ID2, ID1
                TypeInd = P.GetTypeInd(SType1, SType2)
                this = (P, TypeInd)
                #add potential to running lists for right bond types
                if not this in PIndex[ID1][ID2]:
                    PIndex[ID1][ID2].append(this)
    #add potentials and remember unique ones as you go
    UniquePairPotentials = []
    for ID1 in range(NTypes):
        for ID2 in range(0, ID1+1):
            PList = PIndex[ID1][ID2]
            if PList:
                if not PList in UniquePairPotentials:
                    UniquePairPotentials.append(PList)
                Name = PName(Sys, PList)
                if HasLocalRho or HasEwald:
                    s += "pair_coeff      %5d %5d %s %s %s\n" % (ID2+1, ID1+1, d['TABLESTYLE'], PairFile, Name)
                else:
                    s += "pair_coeff      %5d %5d %s %s\n" % (ID2+1, ID1+1, PairFile, Name)
            else:
                if HasLocalRho or HasEwald:
                    s += "pair_coeff      %5d %5d %s %s %s\n" % (ID2+1, ID1+1, d['TABLESTYLE'], PairFile, 'NonBondNull')
                else:
                    s += "pair_coeff      %5d %5d %s %s\n" % (ID2+1, ID1+1, PairFile, 'NonBondNull')
    
    #add local density potentials
    if HasLocalRho:
        LD = []
        for (i, P) in enumerate(Sys.ForceField):
            if not P.Type == ptypes.LocalDensityPotential: continue
            CenterIDs = []
            NeighIDs = []
            for (SType1, SType2) in P.Filter.Select(STypes):
                ID1, ID2 = GetID(SType1), GetID(SType2)
                if not ID1 in CenterIDs: CenterIDs.append(ID1)
                if not ID2 in NeighIDs: NeighIDs.append(ID2)
            LD.append((P, CenterIDs, NeighIDs))
        s += "pair_coeff      * * localdensity %s\n" % LocalDensityFile
        MakeLocalDensityFile(LocalDensityFile, LD)
                
    #add bonded potentials and record 
    AllBonds = []
    AllAngles = []
    AllKeys = []
    NeedsShake = False
    NBondTypes = 0
    NAngleTypes = 0
    RBondTypes = []
    RAngleTypes = []
    RBFConst = RBondFConst * Sys.Units.EScale * Sys.Units.Angstrom**(-2)
    
    #check for rigid molecule
    for MType in Sys.World:
        if not any([m.MType == MType for m in Sys.Mol]): continue
        if MType.Rigid:
            NeedsShake = True
            n = len(MType)
            if n > 4:
                raise LammpsError("Can only format LAMMPS shake input for rigid 2, 3, or 4 site molecules.")
            #enumerate rigid bonds to central molecule
            Pos0 = MType.COMPos[0]
            for i in range(1,n):
                Posi = MType.COMPos[i]
                #add the bond
                bl = sim.geom.Length(Pos0 - Posi)
                s += "bond_coeff      %5d  %-11.4e %-11.4e\n" % (NBondTypes+1, RBFConst, bl)
                RBondTypes.append(NBondTypes)
                #find all of the instances in the system
                for (m, Mol) in enumerate(Sys.Mol):
                    if not Mol.MType is MType: continue
                    aindkey = AtomsIndexKey([Mol[0], Mol[i]])
                    AllBonds.append((aindkey, NBondTypes, True)) 
                NBondTypes += 1
            #enumerate rigid angles to central molecule
            Pos0 = MType.COMPos[0]
            for i in range(1,n):
                Posi = MType.COMPos[i]
                for j in range(i+1,n):
                    Posj = MType.COMPos[j]
                    #add the angle
                    ang = sim.geom.Angle(Posi, Pos0, Posj) * 180. / np.pi
                    Name = "RigidAng%d" % (NAngleTypes + 1)
                    s += "angle_coeff     %5d  %s %s\n" % (NAngleTypes+1, AngleFile, Name)
                    RAngleTypes.append((NAngleTypes, ang))
                    #find all of the instances in the system
                    for (m, Mol) in enumerate(Sys.Mol):
                        if not Mol.MType is MType: continue
                        aindkey = AtomsIndexKey([Mol[i], Mol[0], Mol[j]])
                        AllAngles.append((aindkey, NAngleTypes, True))   
                    NAngleTypes = NAngleTypes + 1
        elif MType.HasRigid:
            raise LammpsError("Can only format LAMMPS shake input for rigid 2, 3, or 4 site molecules.")
            
    #make file for pair table potentials 
    s2 = ""
    for PList in UniquePairPotentials:
        s2 += PotentialTablePair(Sys, PList, GlobalCut)
        s2 += "\n\n"
    #add a null pair potential
    s2 += PotentialTableNullPair(GlobalCut)
    s2 += "\n\n"
    #write to file
    file(PairFile, "w").write(s2) 
    
    #do angle potentials
    s2 = ""
    for P in Sys.ForceField:
        if not P.Type == ptypes.AnglePotential: continue
        s2 += PotentialTableAngle(Sys, P)
        s2 += "\n\n"
    #add rigid angles
    for (Ind, Theta0) in RAngleTypes:
        Name = "RigidAng%d" % (Ind+1)
        s2 += PotentialTableAngleRigid(Sys, Name, Theta0)
        s2 += "\n\n"
    #write to file
    if len(s2):
        file(AngleFile, "w").write(s2)

    #make potentials for dihedrals
    s2 = ""
    for P in Sys.ForceField:
        if not P.Type == ptypes.TorsionPotential: continue
        s2 += PotentialTableDihedral(Sys, P)
        s2 += "\n\n"
    #write to file
    if len(s2):
        file(DihedralFile, 'w').write(s2)        
                
    #now do harmonic bonds or spline bonds
    BondedPotentials = []
    for P in Sys.ForceField:
        if isinstance(P, sim.potential.Ewald): continue
        if P.Type == ptypes.PairPotential and P.Filter.Bonded:
            BondedPotentials.append(P)
            #select atoms that work for this potential
            FoundBond = False
            for Atoms in P.Filter.Select(Sys.Atom):
                aindkey = AtomsIndexKey(Atoms)
                if aindkey in AllKeys:
                    raise LammpsError("Found multiple bond potentials for atoms " + str(aindkey))
                FoundBond = True
                AllBonds.append((aindkey, NBondTypes, False))
                AllKeys.append(aindkey)
            if FoundBond:
                if BondSpline:
                    Name = PName(Sys, P)  
                    s += "bond_coeff      %5d  %s %s\n" % (NBondTypes+1, BondFile, Name)
                else:
                    s += "bond_coeff      %5d  %-11.4e %-11.4e\n" % (NBondTypes+1, P.FConst[0], P.Dist0[0])                
                NBondTypes += 1
    AllBonds.sort()               
    del AllKeys

    #make file for bond table potentials
    if BondSpline:
        s2 = ""
        for P in BondedPotentials:
            BCut = P.Cut
            s2 += PotentialTableBond(Sys, P, BCut)
            s2 += "\n\n"
        #write to file
        file(BondFile, "w").write(s2)

    #add angle potentials and record 
    AllKeys = []
    for P in Sys.ForceField:
        if P.Type == ptypes.AnglePotential:
            #select atoms that work for this potential
            FoundAngle = False
            for Atoms in P.Filter.Select(Sys.Atom):
                FoundAngle = True
                aindkey = AtomsIndexKey(Atoms)
                if aindkey in AllKeys:
                    raise LammpsError("Found multiple angle potentials for atoms " + str(aindkey))
                AllAngles.append((aindkey, NAngleTypes, False))
                AllKeys.append(aindkey)
                Name = PName(Sys, P)
            if FoundAngle:
                s += "angle_coeff     %5d  %s %s\n" % (NAngleTypes+1, AngleFile, Name)
                NAngleTypes += 1
    AllAngles.sort()
    del AllKeys
    
    #add dihedral potentials and record 
    AllDihedrals = []
    AllKeys = []
    NDihedralTypes = 0
    for P in Sys.ForceField:
        if P.Type == ptypes.TorsionPotential:
            #select atoms that work for this potential
            FoundDihedral = False
            for Atoms in P.Filter.Select(Sys.Atom):
                FoundDihedral = True
                aindkey = AtomsIndexKey(Atoms)
                if aindkey in AllKeys:
                    raise LammpsError("Found multiple dihedral potentials for atoms " + str(aindkey))
                AllDihedrals.append((aindkey, NDihedralTypes))
                AllKeys.append(aindkey)
                Name = PName(Sys, P)
            if FoundDihedral:
                s += "dihedral_coeff  %5d  %s %s\n" % (NDihedralTypes+1, DihedralFile, Name)
                NDihedralTypes += 1
    AllDihedrals.sort()
    del AllKeys

    #add the special bonds coefficients
    s += "\n"
    s += 'special_bonds     lj %3.1f %3.1f %3.1f' % (LJCoef12, LJCoef13, LJCoef14)
    if HasEwald:
        s += ' coul %3.1f %3.1f %3.1f' % (EWCoef12, EWCoef13, EWCoef14)
    s += "\n"
    
    #add the shake command
    if NeedsShake:
        RigidCmd = "fix shakefix all shake %11.4e %d 0 b " % (ShakeTol, ShakeIter)
        RigidCmd = RigidCmd + " ".join(["%d" % (i+1) for i in RBondTypes]) 
        if len(RAngleTypes):
            RigidCmd += " a " + " ".join(["%d" % (i+1) for (i,ang) in RAngleTypes])
    else:
        RigidCmd = ""
        
    #add the skin distance
    if not SkinDistance is None:
        s += "neighbor %11.3e bin\n" % SkinDistance
    #add the neighbor commands
    if not NeighDelay is None:
        s += "neigh_modify delay %d\n" % NeighDelay
    if not NeighEvery is None:
        s += "neigh_modify every %d\n" % NeighEvery
    if not NeighCheck is None:
        s += "neigh_modify check %s\n" % NeighCheck       
    if not NeighPage is None: 
        s += "neigh_modify page %d\n" % NeighPage
    if not NeighOne is None:
        s += "neigh_modify one  %d\n" % NeighOne
     
    #add external potentials -- need to manually add functional form of potential
    if HasExtField:
        ExtFields = []

        s += "\n"
        ImplementedExtFields = ["external_sinusoid"]
        for (i, P) in enumerate(Sys.ForceField):
            if not P.Type == ptypes.FieldPotential: continue
            if P.Names[0] not in ImplementedExtFields: continue
           
            ExtFields.append(P)
            s += P.LammpsStr() + "\n"        

    #add any lammps commands
    if len(LammpsCommands):
        LammpsCommands = LammpsCommands.replace("[RIGIDFIX]", RigidCmd)
        s += "\n\n" + LammpsCommands % d

    #add the atom names
    s = s.replace("ATOMNAMES", AtomNames)

    #write out the input file
    file(InFile, "w").write(s)    
    
    #make a header for the data file    
    t = (len(Sys.Atom), len(AllBonds), len(AllAngles), len(AllDihedrals), 0,
         NTypes, NBondTypes, NAngleTypes, NDihedralTypes, 0)
    s1 = """LAMMPS Description

%12d  atoms
%12d  bonds
%12d  angles
%12d  dihedrals
%12d  impropers

%12d  atom types
%12d  bond types
%12d  angle types
%12d  dihedral types
%12d  improper types
""" % t

    #get the atomic positions
    Pos = Sys.Pos.copy()
    #get the minimum image translations
    nImage = sim.geom.MinimageN(Pos, Sys.BoxL)
    #reset the positions inside of the box
    Pos = sim.geom.Minimage(Pos, Sys.BoxL)

    #add the box dimensions
    for (i, ax) in enumerate(['x', 'y', 'z']):      
        if Sys.BoxL[i] > 0:
            lo, hi = -0.5 * Sys.BoxL[i], 0.5 * Sys.BoxL[i]
        else:
            #center along this dimension
            Pos[:,i] = Pos[:,i] - Pos[:,i].mean()
            #now find the min and max position
            maxpos = np.max(Pos[:,i])
            minpos = np.min(Pos[:,i])
            #find the maximum atom extent and then double for safety 
            #(LAMMPS will shrink-wrap anyways)
            h = max(abs(maxpos), abs(minpos)) * 2.
            lo, hi = -h, h
        s1 += "    %11.4e %11.4e %slo %shi\n" % (lo, hi, ax, ax)
    s1 += "\n"
    
    #add the masses
    s1 += 'Masses\n\n'
    for Type1 in Types:
        s1 += '%3d %11.4f\n' % (GetID(Type1)+1, Type1.Mass)
    s1 += "\n"

    #add the atoms
    #format is atom-ID molecule-ID atom-type x y z for unchaarged
    # and atom-ID molecule-ID atom-type charge x y z for charged
    s1 += 'Atoms\n\n'
    for (i, Atom) in enumerate(Sys.Atom):
        aPos = Pos[i]
        aCharge = Atom.Charge
        aN = -nImage[i]
        if ResetPosFirst:
            aPos = sim.geom.Minimage(aPos, Sys.BoxL)
        if HasEwald:
            s1 += "%10d %5d %5d %10.5f %12.5e %12.5e %12.5e %d %d %d\n" % (Atom.Ind+1, Atom.Parent.MInd+1, GetID(Atom)+1,
                                                                           aCharge, aPos[0], aPos[1], aPos[2],
                                                                           aN[0], aN[1], aN[2])            
        else:
            s1 += "%10d %5d %5d %12.5e %12.5e %12.5e %d %d %d\n" % (Atom.Ind+1, Atom.Parent.MInd+1, GetID(Atom)+1,
                                                                    aPos[0], aPos[1], aPos[2],
                                                                    aN[0], aN[1], aN[2])
    s1 += "\n"
    
    #add the bonds
    #format is ID type atom1 atom2
    if len(AllBonds):
        s1 += '\nBonds\n\n'
        for (ID, ((AInd1, AInd2), TypeInd, Rigid)) in enumerate(AllBonds):
            s1 += "%10d %10d %10d %10d\n" % (ID+1, TypeInd+1, AInd1+1, AInd2+1)
        s1 += "\n"
    
    #add the angles
    #format is ID type atom1 atom2 atom3
    if len(AllAngles):
        s1 += '\nAngles\n\n'
        for (ID, ((AInd1, AInd2, AInd3), TypeInd, Rigid)) in enumerate(AllAngles):
            s1 += "%10d %10d %10d %10d %10d\n" % (ID+1, TypeInd+1, AInd1+1, AInd2+1, AInd3+1)
        s1 += "\n"
    
    #add the dihedrals
    #format is ID type atom1 atom2 atom3 atom3
    if len(AllDihedrals):
        s1 += '\nDihedrals\n\n'
        for (ID, ((AInd1, AInd2, AInd3, AInd4), TypeInd)) in enumerate(AllDihedrals):
            s1 += "%10d %10d %10d %10d %10d %10d\n" % (ID+1, TypeInd+1, AInd1+1, AInd2+1, AInd3+1, AInd4+1)
        s1 += "\n"
    
    #add the atom names
    s1 = s1.replace("ATOMNAMES", AtomNames)

    #write out the data file
    file(DataFile, "w").write(s1)
    
    #return the filenames, first is the input file
    LammpsFiles = [InFile, DataFile, PairFile]
    if HasAngle:
        LammpsFiles.append(AngleFile)
    if HasDihedral:
        LammpsFiles.append(DihedralFile)
    if HasLocalRho:
        LammpsFiles.append(LocalDensityFile)
    return LammpsFiles


def MakeLammpsMD(Sys, NStepsMin = 10000, NStepsEquil = 100000, NStepsProd = 5000000, 
                 WriteFreq = 1000, Prefix = "", TrajFile = "trj.lammpstrj",
                 WriteRestart = True, RestartFile = None, 
                 ThermoOutput = "pe temp press", ConstPressure = False, UseLangevin = False,
                 LammpsCommandsBefore = "", LammpsCommandsAfter = "", OutputDCD = False,
                 CalcPress = False, Nevery=1, Nrepeat=1, Nfreq=1, ThermoStyle = "line",*args, **kwargs):
    """Makes a LAMMPS MD trajectory, with default settings.
Returns InFile, DataFile, Pairfile, AngleFile, DihedralFile, TrajFile"""
    import sim.units
    
    #add a prefix
    TrajFile = Prefix + TrajFile 
    PressFile = Sys.Name+'_PressData/'+'Pressure.dat'
    RgFile    = Sys.Name+'_RgData'
    Restart1 = Sys.Name+ '_1.restart'
    Restart2 = Sys.Name+ '_2.restart'
    #start input file
    s = ""
    if len(LammpsCommandsBefore):
        s += LammpsCommandsBefore + "\n\n"

    #choose the timestep 
    if Sys.Units == sim.units.DimensionlessUnits:
        TimeStep = Sys.Int.Methods.VVIntegrate.TimeStep
    elif Sys.Units == sim.units.MKSUnits:
        TimeStep = Sys.Int.Methods.VVIntegrate.TimeStep
    elif Sys.Units == sim.units.AtomicUnits:
        TimeStep = Sys.Int.Methods.VVIntegrate.TimeStep * 1000. / 20.5  #convert to fs
    else:
        raise LammpsError("Don't recognize type of units in system.")
    #check integrator
    if not Sys.Int.Method.Barostat == Sys.Int.Method.BarostatNone:
        ConstPressure = True
    if Sys.Int.Method.Thermostat ==  Sys.Int.Method.ThermostatLangevin:
        UseLangevin = True                    
    #make a dictionary for filling in template
    RanSeed = 100000 * np.random.rand()
    d = {"TIMESTEP" : TimeStep,
         "TEMP" : Sys.TempSet,
         "PRES" : Sys.PresSet,
         "LANGEVINDAMP" : 1. / Sys.Int.Methods.VVIntegrate.LangevinGamma,
         "TEMPDAMP" : 100. * TimeStep,
         "PRESDAMP" : PDampFactor * TimeStep,
         "WRITEFREQ" : WriteFreq,
         "NSTEPSMIN" : NStepsMin,
         "NSTEPSEQUIL" : NStepsEquil,
         "NSTEPSPROD" : NStepsProd,
         "TRAJFILE" : TrajFile,
         "RANSEED" : RanSeed,
         "THERMOOUTPUT" : ThermoOutput,
         "PRESSFILE" : PressFile,
         "NUMMOL" : Sys.NMol,
         "NEVERY" : Nevery,
         "NREPEAT" : Nrepeat,
         "Nfreq" : Nfreq,
         "RGFILE" : RgFile,
         "RESTART1": Restart1,
         "RESTART2": Restart2
         }
         
    #add to input string
    if not ThermoStyle == 'line':
        s += """
#Thermostat & time integration
timestep        %(TIMESTEP)-11.4e
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   norm no format float %%%%14.7e

#minimization
minimize        1.e-4 0.0 %(NSTEPSMIN)d %(NSTEPSMIN)d
""" % d
    else:
        s += """
#Thermostat & time integration
timestep        %(TIMESTEP)-11.4e
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   line multi norm no format float %%%%14.7e

#minimization
minimize        1.e-4 0.0 %(NSTEPSMIN)d %(NSTEPSMIN)d
""" % d

    #check if NVT or NPT
    if ConstPressure:
        if UseLangevin:
            raise ValueError("Cannot use Langevin thermostat with NPT in LAMMPS export.")
        s += """
#Set up integration parameters
fix             nptstat all npt temp %(TEMP)11.4e %(TEMP)11.4e %(TEMPDAMP)11.4e iso %(PRES)11.4e %(PRES)11.4e %(PRESDAMP)11.4e
[RIGIDFIX]
""" % d
    elif UseLangevin:
        s += """
#set up integration parameters
fix             timeintegration all nve
fix             thermostat all langevin %(TEMP)11.4e %(TEMP)11.4e %(LANGEVINDAMP)11.4e %(RANSEED)d
[RIGIDFIX]
""" % d  
    else:       
        s += """
#set up integration parameters        
fix             nvtstat all nvt temp %(TEMP)11.4e %(TEMP)11.4e %(TEMPDAMP)11.4e
[RIGIDFIX]
""" %d 

    #recentering if not periodic
    if np.any(Sys.BoxL <= 0):
        s += "fix             recenterfix all recenter 0.0 0.0 0.0 units box\n"
        
        #run equilibration and production
    if not ThermoStyle == 'line':
        s += """
#restart thermo
reset_timestep  0
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   norm no format float %%%%14.7e
variable PRESS equal press
variable STEP equal step

#run equilibration
run             %(NSTEPSEQUIL)d

#restart thermo
reset_timestep  0
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   norm no format float %%%%14.7e

""" %d
    else:
        s += """
#restart thermo
reset_timestep  0
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   line multi norm no format float %%%%14.7e
variable PRESS equal press
variable STEP equal step

#run equilibration
run             %(NSTEPSEQUIL)d

#restart thermo
reset_timestep  0
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   line multi norm no format float %%%%14.7e
""" %d 

    if CalcPress:
        s += """
fix print_press all print %(WRITEFREQ)d "${STEP} ${PRESS}" file %(PRESSFILE)s screen no title "# STEP PRESS"
""" %d

    if OutputDCD:
        s += """
#output .DCD file instead of lammpstrj
dump            myDump all dcd %(WRITEFREQ)d %(TRAJFILE)s 
""" % d

    else:
        s += """
#setup trajectory output
dump            myDump all custom %(WRITEFREQ)d %(TRAJFILE)s id type x y z mol
dump_modify     myDump element ATOMNAMES
dump_modify     myDump sort id
dump_modify     myDump format line  "%%%%d %%%%d %%%%14.7e %%%%14.7e %%%%14.7e %%%%d"
""" % d

    if WriteRestart:
        s += """
restart 1000 %(RESTART1)s %(RESTART2)s""" % d

    s += """
# PRODUCTION RUNS
#run production
run             %(NSTEPSPROD)d
""" % d

    #recentering if not periodic
    if np.any(Sys.BoxL <= 0):
        s += "unfix           recenterfix\n"
    
    if len(LammpsCommandsAfter):
        s += "\n\n" + LammpsCommandsAfter + "\n"
        
    #now make the lammps files
    LammpsFiles = MakeLammps(Sys, *args, Prefix = Prefix, LammpsCommands = s, RestartFile = RestartFile, **kwargs)
    #return filenames
    return LammpsFiles, TrajFile
    

def CheckLammpsOutput(LogFile):
    """Parses LAMMPS output to check for any problems."""
    import os
    if not os.path.isfile(LogFile):
        raise LammpsError("Could not find LAMMPS log file %s" % LogFile)
    f = file(LogFile)
    s = f.readline()
    NeighPct = None    
    while len(s):
        if s.startswith("Neigh time (%) ="):
            try:
                NeighPct = float(s.split()[-1][1:-1])
            except StandardError:
                pass
        elif s.startswith("Neigh   |"):
                NeighPct = float(s.split("|")[-1])
        s = f.readline()
    f.close()
    if not NeighPct is None:
        if WARNINGS and NeighPct > 25:
            print "WARNING: LAMMPS spent over %.1f%% of time making neighbor lists." % NeighPct
            print "         May want to increase SkinDistance."
            
   
def RunLammps(InFile, Prefix = "", LogFile = "lammps.log", Verbose = False,
              CheckOutput = True):
    """Runs lammps.  Returns LogFile, ReturnCode."""
    import subprocess, os
    import time
    #add prefix
    LogFile = Prefix + LogFile
    args = [LammpsExec, "-in", InFile, "-log", LogFile]
    if Verbose:
        print "Running LAMMPS with command:"
        cmd = " ".join(args)
        print("{}".format(cmd))
    t1 = time.time()
    #send screen output to null to prevent memory overflow
    DEVNULL = open(os.devnull, "w")
    p = subprocess.Popen(cmd, shell = True, 
                         stdout = DEVNULL,
                         stderr = subprocess.PIPE)
    retout, reterr = p.communicate()
    returncode = p.returncode
    if Verbose:
        t2 = time.time()
        print "LAMMPS finished in %d seconds with return code %d." % (t2-t1, returncode)
    if not returncode == 0:
        #there was an error
        print "="*20 + " LAMMPS OUTPUT " + "="*20
        if not retout is None: print retout
        if not reterr is None: print reterr
        print "="*20 + " END LAMMPS OUTPUT " + "="*20
        print "Input file: %s" % InFile
        raise LammpsError("Nonzero return code found when running LAMMPS.")
    if CheckOutput:
        CheckLammpsOutput(LogFile)
    return LogFile, returncode




def MakeLammpsTraj(Sys, ReturnTraj = True, DelTempFiles = True, Prefix = "", Verbose = False, *args, **kwargs):
    """Returns a LAMMPS trajectory object.  Takes same arguments as MakeLammpsMD.
Returns Traj, TrajFile."""
    import sim.traj
    import os
    LammpsFiles, TrajFile = MakeLammpsMD(Sys, *args, Prefix = Prefix, **kwargs)
    InFile = LammpsFiles[0]
    LogFile, returncode = RunLammps(InFile, Prefix = Prefix, Verbose = Verbose)
    if ReturnTraj:
        Traj = sim.traj.Lammps(TrajFile, LogFile = LogFile, LogFileToken = "#run production")
    else:
        Traj = 'N/A'
    if DelTempFiles:
        for fn in LammpsFiles + [LogFile]:
            if os.path.isfile(fn):
                os.remove(fn)
    return Traj, TrajFile

def MakeRerunLammpsTraj(Sys, ElecTrajFile, ReturnTraj = True, DelTempFiles = True, Prefix = "", Verbose = False, *args, **kwargs):
    """Reprocess trajectory using forcefield of Sys
Returns a LAMMPS trajectory object.  Takes same arguments as MakeLammpsMD.
Returns Traj, TrajFile.
ElecTrajFile: trajectory to reprocess with forcefield from Sys"""
    import sim.traj
    import os
    LammpsFiles = MakeRerunLammpsMD(Sys, ElecTrajFile,*args, Prefix = Prefix, **kwargs)
    InFile = LammpsFiles[0]
    print('\nReprocess trajectory')
    LogFile, returncode = RunLammps(InFile, Prefix = Prefix, Verbose = Verbose)
    if ReturnTraj:
        Traj = sim.traj.Lammps(ElecTrajFile, LogFile = LogFile, LogFileToken = "#run production")
    else:
        Traj = 'N/A'
    if DelTempFiles:
        for fn in LammpsFiles + [LogFile]:
            if os.path.isfile(fn):
                os.remove(fn)
    return Traj


def MakeRerunLammpsMD(Sys, ElecTrajFile, NStepsMin = 10000, NStepsEquil = 100000, NStepsProd = 5000000, 
                 WriteFreq = 1000, Prefix = "", TrajFile = "trj.lammpstrj",
                 WriteRestart = True, RestartFile = None, 
                 ThermoOutput = "pe temp press", ConstPressure = False, UseLangevin = False,
                 LammpsCommandsBefore = "", LammpsCommandsAfter = "", OutputDCD = False,
                 CalcPress = False, Nevery=1, Nrepeat=1, Nfreq=1, ThermoStyle = "line",*args, **kwargs):
    """Makes a LAMMPS MD trajectory, with default settings.
Returns InFile, DataFile, Pairfile, AngleFile, DihedralFile, TrajFile"""
    import sim.units
    
    #start input file
    s = ""
    if len(LammpsCommandsBefore):
        s += LammpsCommandsBefore + "\n\n"

    #choose the timestep 
    if Sys.Units == sim.units.DimensionlessUnits:
        TimeStep = Sys.Int.Methods.VVIntegrate.TimeStep
    elif Sys.Units == sim.units.MKSUnits:
        TimeStep = Sys.Int.Methods.VVIntegrate.TimeStep
    elif Sys.Units == sim.units.AtomicUnits:
        TimeStep = Sys.Int.Methods.VVIntegrate.TimeStep * 1000. / 20.5  #convert to fs
    else:
        raise LammpsError("Don't recognize type of units in system.")
    #check integrator
    if not Sys.Int.Method.Barostat == Sys.Int.Method.BarostatNone:
        ConstPressure = True
    if Sys.Int.Method.Thermostat ==  Sys.Int.Method.ThermostatLangevin:
        UseLangevin = True                    
    #make a dictionary for filling in template
    RanSeed = 100000 * np.random.rand()
    d = {"TIMESTEP" : TimeStep,
         "TEMP" : Sys.TempSet,
         "PRES" : Sys.PresSet,
         "LANGEVINDAMP" : 1. / Sys.Int.Methods.VVIntegrate.LangevinGamma,
         "TEMPDAMP" : 100. * TimeStep,
         "PRESDAMP" : PDampFactor * TimeStep,
         "WRITEFREQ" : 1, #always 1 for rerun so that report energy from every frame
         "NSTEPSMIN" : NStepsMin,
         "NSTEPSEQUIL" : NStepsEquil,
         "NSTEPSPROD" : NStepsProd,
         "TRAJFILE" : TrajFile,
         "RANSEED" : RanSeed,
         "THERMOOUTPUT" : ThermoOutput,
         "NUMMOL" : Sys.NMol,
         "NEVERY" : Nevery,
         "NREPEAT" : Nrepeat,
         "Nfreq" : Nfreq,
         "ELECTRAJFILE": ElecTrajFile,
         }
         
    #add to input string
    if not ThermoStyle == 'line':
        s += """
#Thermostat & time integration
timestep        %(TIMESTEP)-11.4e
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   norm no format float %%%%14.7e
""" % d
    else:
        s += """
#Thermostat & time integration
timestep        %(TIMESTEP)-11.4e
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   line multi norm no format float %%%%14.7e
""" % d

    #check if NVT or NPT
    if ConstPressure:
        if UseLangevin:
            raise ValueError("Cannot use Langevin thermostat with NPT in LAMMPS export.")
        s += """
#Set up integration parameters
fix             nptstat all npt temp %(TEMP)11.4e %(TEMP)11.4e %(TEMPDAMP)11.4e iso %(PRES)11.4e %(PRES)11.4e %(PRESDAMP)11.4e
[RIGIDFIX]
""" % d
    elif UseLangevin:
        s += """
#set up integration parameters
fix             timeintegration all nve
fix             thermostat all langevin %(TEMP)11.4e %(TEMP)11.4e %(LANGEVINDAMP)11.4e %(RANSEED)d
[RIGIDFIX]
""" % d  
    else:       
        s += """
#set up integration parameters        
fix             nvtstat all nvt temp %(TEMP)11.4e %(TEMP)11.4e %(TEMPDAMP)11.4e
[RIGIDFIX]
""" %d 

    #recentering if not periodic
    if np.any(Sys.BoxL <= 0):
        s += "fix             recenterfix all recenter 0.0 0.0 0.0 units box\n"
        
    if not ThermoStyle == 'line':
        s += """
#restart thermo
reset_timestep  0
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   norm no format float %%%%14.7e
""" %d
    else:
        s += """
#restart thermo
reset_timestep  0
thermo          %(WRITEFREQ)d
thermo_style    custom %(THERMOOUTPUT)s
thermo_modify   line multi norm no format float %%%%14.7e
""" %d 


    s += """
# PRODUCTION RUNS
#run production
rerun %(ELECTRAJFILE)s dump x y z 
""" %d

    #recentering if not periodic
    if np.any(Sys.BoxL <= 0):
        s += "unfix           recenterfix\n"
    
    if len(LammpsCommandsAfter):
        s += "\n\n" + LammpsCommandsAfter + "\n"
        
    #now make the lammps files
    LammpsFiles = MakeLammps(Sys, *args, Prefix = Prefix, LammpsCommands = s, RestartFile = RestartFile, **kwargs)
    #return filenames
    return LammpsFiles

