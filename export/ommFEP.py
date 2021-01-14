#/usr/bin/env python

### Translates Sys objects to LAMMPS input files, in SIM suite.
### coded by KS, MSS
### Overall code flow:
#   1) optimizetraj calls MakeOpenMMTraj(...)
#   2) MakeOpenMMTraj manages files and calls RunOpenMM(...)
#   3) RunOpenMM(...) contains the minimization, equilibration, production run, and calls CreateOMMSimulation(),CreateOMMSys()
#   4) CreateOMMSimulation(...) creates simulation object (i.e. integrator and barostat)
#   5) CreateOMMSys(...) does the heavy work of creating topology, system, setting up force fields

### Wishlist to be implemented:
#   1) electrostatics
#   2) other bonding potentials
#   3) bonding exclusions
#   4) angle constraints

import os

import numpy as np

import sim.potential.base.potentialtypes as ptypes
import sim.spline as spline

#========================================
#Platform settings
#========================================
device = -1 #-1 to let omm automatically choose platform type, device, and precision
platformName = 'CUDA'
platformPrecision = 'mixed' #activated only if device != -1

#========================================
#Some simulation settings
#========================================
constraintTolerance = 1e-6

#========================================
#Debug options
#========================================
WARNINGS = True
VVerbose = False #Very Verbose


#========================================
#Other settings
#========================================
UseTabulated = True
NPairPotentialKnots = 251 #note every knot point has 4 spline coefficients
NPairPotentialBins = 1001 
NAnglePotentialBins = 500
NDihedralPotentialBins = 500
NLocalDensityPotentialBins = 500
NOtherPotentialBins = 500

#pair spline inner cutoff in length scale
InnerCutoff = 0.001

#max pair energy scale in units of kB T
MaxPairEnekBT = 20.


#========================================
#Legacy settings from lammps, unused, possibly delete
#========================================
#neighbor list settings
NeighDelay = None
NeighEvery = None
NeighCheck = None
SkinDistance = None

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



# === make an error class ===
class OMMError(Exception):
  def __init__(self, Msg):
    self.Msg = Msg
  def __str__(self):
    return str(self.Msg)
 
#=== Helper functions adapted from Lammps Export ===
def TabulatePair(Sys,PList,Cut):
    """Takes a list of pair potentials, and returns tabulated values, using global discretization settings.
    OpenMM's ContinuousFunctions take evenly spaced points, so InnerCutoff setting is only approximate"""
    #HERE RAISE FLAG IF NO CUTOFF
    if Cut is None:
        print("Automatically choosing maximum cutoff")
        xmax = np.max([P.Arg.ReportMax.max() for (P, TypeInd) in PList])
    else:
        xmax = Cut
   
    global NPairPotentialBins
    N = NPairPotentialBins
    dxguess = xmax/N
    if dxguess <= InnerCutoff: 
        xmin = dxguess
    else:
        print("Caution: Increasing NPairPotential to approach InnerCutoff! Assuming this is Global, i.e. all potentials using same Cutoff")
        N = np.ceil(Cut/InnerCutoff)
        NPairPotentialBins = N
        xmin = xmax/N
    xvals = np.linspace(xmin,xmax,N)


    uvals = np.zeros_like(xvals)
    duvals = np.zeros_like(xvals)
    
    for (P, TypeInd) in PList:
        P.SetTypeInd(TypeInd)
        uvals += np.array([P.Val(x) for x in xvals])
        duvals += np.array([-P.DVal(x) for x in xvals]) #actually the negative of the derivative, i.e. force
    if not np.isfinite(uvals[0]): uvals[0] = 1.e300

    #get the maximum pair energy
    #MaxPairEne = Sys.Units.kB * Sys.TempSet * MaxPairEnekBT
    MaxPairEne = MaxPairEnekBT #right now default is dimensionless units

    #indices where energy is greater
    ind = np.where(uvals > MaxPairEne)[0]
    if len(ind):
        #find the first index where energy is valid
        i = ind[-1] + 1
        #do a linear extrapolation in the hard core region
        uvals[:i] = (xvals[i] - xvals[:i]) * duvals[i] + uvals[i]
        duvals[:i] = duvals[i]

    #insert value at x=0 so that still valid if r<InnerCutoff
    xvals = np.insert(xvals,0,0.0)
    uvals = np.insert(uvals,0,uvals[0])
    duvals = np.insert(duvals,0,0.0)

    return xvals,uvals,xmax

def GetKnots(rcut,nknots,xfit,yvals):
    sp = spline.Spline(rcut,np.zeros(nknots))
    sp.fitCoeff(xfit,yvals)
    knots = sp.knots
    spCoeffs = sp.SPCoeff
    return knots,spCoeffs

# ===================== #
# OpenMM Implementation #
# ===================== #
try:
    from simtk import openmm, unit
    #from simtk.unit import *
    from simtk.openmm import app
    import numpy as np
    import simtk.openmm.app.topology as topology
    import mdtraj

    import sys, time, os
    cwd=os.path.dirname( os.path.realpath(__file__) )
    simpath="/".join(cwd.split("/")[:-3])
    if sys.path[1] != simpath:
        sys.path.insert(1,simpath)
    #print(sys.path)

    import sim
    import sim.potential.base.potentialtypes as ptypes



    #========================================
    ###DEFINE ENERGY, LENGTH & MASS SCALES###
    #Try to use OMM's built-in unit definitions to the extent possible#
    #========================================
    epsilon  = 1.0 * unit.kilojoules_per_mole     #kJ/mole
    sigma = 1.0 * unit.nanometer                #nm
    tau = 1.0*unit.picoseconds                  #ps
    mass = tau**2/sigma**2*epsilon              #dalton = g/mole

    #N_av = 6.022140857*10**23 /unit.mole
    #kb = 1.380649*10**(-23)*unit.joules/unit.kelvin* N_av #joules/kelvin/mol
    N_av = unit.constants.AVOGADRO_CONSTANT_NA #6.02214179e23 / mole
    kB = unit.constants.BOLTZMANN_CONSTANT_kB #1.3806504e-23 * joule / kelvin
    Tref = epsilon/kB/N_av

    #from openmmtools.constants import ONE_4PI_EPS0 
    ONE_4PI_EPS0 = 138.935456       #in OMM units, 1/4pi*eps0, [=] kT/mole*nm/e^2
    qFactor = ONE_4PI_EPS0**-0.5    #s.t. electrostatic energy (1/4pi eps0)*q^2*qF^2/r = 1kB*Tref*Nav/mole = 1kJ/mole = 1epsilon

except:
    print("openmm or mdtraj did not load; these features will not work!")


def CreateOMMSys(Sys,Verbose=False):
    print(Sys.Units.Name)
    if Sys.Units.Name != "DimensionlessUnits": raise ValueError( "Danger! OMM export currently only for dimensionless units, but {} detected".format(Sys.Units.Name) )
    "takes in Sim 'Sys' object and returns an OpenMM System and Topology"
    print("=== OMM export currently uses LJ (Dimensionless) Units ===")
    print("mass: {} dalton".format(mass.value_in_unit(unit.dalton)))
    print("epsilon: {}".format(epsilon))
    print("tau: {}".format(tau))

    # --- box size ---
    Lx, Ly, Lz = Sys.BoxL
   
    # --- Other Options ---
    # TODO: nonbonded method, electrostatics, bond constraints
    #nonbondedMethod = openmm.CustomNonbondedForce.CutoffPeriodic
    #ewaldErrorTolerance = 0.0001
    #constraints = None
    #nonbondedCutoff = reduced_nonbondedCutoff*sigma

    #========================================
    # Create a system and add particles to it
    #========================================
    print("=== Creating OMM System ===")
    system = openmm.System()

    # Set the periodic box vectors:
    box_edge = [Lx,Ly,Lz]
    box_vectors = np.diag(box_edge) * sigma
    system.setDefaultPeriodicBoxVectors(*box_vectors)
#    print('Box vectors:\n{}'.format*(system.getPeriodicBoxVectors()))
    #==================
    ##CREATE TOPOLOGY##
    #==================
    #Topology consists of a set of Chains 
    #Each Chain contains a set of Residues, 
    #and each Residue contains a set of Atoms.
    #We take the topology from the `sim System` object. Currently we do 1 residue per chain, and residue is effectively the Molecule class in sim.
    print("\n--- Creating Topology ---")
    # --- first get atom types and list ---
    simAtomTypes = Sys.World.AtomTypes
    elements = {}
    atomTypeIndex = {}
    atomNameMap = []
    simAtomTypeMap = {}
    for ia,a in enumerate(simAtomTypes):
        newsymbol = 'Z{}'.format(ia)
        if newsymbol not in app.element.Element._elements_by_symbol:
            elements[a.Name]=app.element.Element(200+ia, a.Name, newsymbol, mass)
        else:
            elements[a.Name]=app.element.Element.getBySymbol(newsymbol)

        atomTypeIndex[a.Name] = ia
        atomNameMap.append(a.Name)
        simAtomTypeMap[a.Name] = a
        #Check that OMM atom index corresponds to that in sim:
        ia == a.AID

    # --- next get the molecules and bondList ---
    residues = []
    moleculeAtomList = {}
    moleculeBondList = {}
    for im,m in enumerate(Sys.World):
        residues.append(m.Name)
        moleculeAtomList[m.Name] = [a.Name for a in m] 
        moleculeBondList[m.Name] = [ [b.SType1.AInd,b.SType2.AInd] for b in m.Bonds ]

    # --- aggregate stuff, and add atoms to omm system ---
    # Particles are added one at a time
    # Their indices in the System will correspond with their indices in the Force objects we will add later
    atomList = [a.Name for a in Sys.Atom]
    molList = [m.Name for a in Sys.Mol]
    for a in Sys.Atom:
        system.addParticle(mass)
    print("Total number of paricles in system: {}".format(system.getNumParticles()))
    
    # --- the actual work of creating the topology ---
    top = app.topology.Topology()
    mdtrajtop = app.topology.Topology() #so that later can make molecules whole
    constrainedBonds = False
    constraintlengths = []
    for im,mol in enumerate(Sys.Mol):
        chain = top.addChain() #Create new chain for each molecule
        res = top.addResidue(mol.Name,chain)
        mdt_chain = mdtrajtop.addChain() #Create new chain for each molecule
        mdt_res = mdtrajtop.addResidue(mol.Name,mdt_chain)
        
        # add the atoms
        atomsInThisRes = []
        mdt_atomsInThisRes = []
        for atomInd,a in enumerate(mol):
            el = elements[a.Name]
            if atomInd > 0:
                previousAtom = atom
            atom = top.addAtom( a.Name, el, res )
            mdt_atom = mdtrajtop.addAtom( a.Name, mdtraj.element.Element.getByAtomicNumber(atomTypeIndex[a.Name]), mdt_res ) #use a dummy element by matching atomic number == cgAtomTypeIndex
            atomsInThisRes.append(atom)
            mdt_atomsInThisRes.append(mdt_atom)

        # add the bonds
        for bond in mol.Bonds:
            a1 = atomsInThisRes[bond.SType1.AInd] #the atomsInThisRes has the current, newly added atoms of this residue
            a2 = atomsInThisRes[bond.SType2.AInd]
            if VVerbose: print("Adding bond ({},{})".format(bond.SType1.AInd,bond.SType2.AInd))
            #top.addBond( atomsInThisRes[bond.SType1.AInd], atomsInThisRes[bond.SType2.AInd] )
            newbond = top.addBond( a1, a2 )
            mdtrajtop.addBond(a1,a2) #don't worry about adding constraint to the mdtraj topology

            if bond.Rigid:
                constrainedBonds = True
                if VVerbose: print("Adding rigid constraint for {},{}".format(a1,a2))
                system.addConstraint( a1.index, a2.index, bond.RLength ) #constraint uses 0-based indexing of atoms
                constraintlengths.append(bond.RLength)
            else:
                constraintlengths.append(0.0)

    print("Total number of constraints: {}".format(system.getNumConstraints()))   

    #====================
    ##CREATE FORCEFIELD##
    #====================
    #Currently only allow for some restricted interactions
    #a) harmonic bonds
    #b) Gaussian repulsion
    #c) external force

    #TODO:
    #1) charges
    #2) spline/tabulated
    #3) special bonds treatment (i.e. all fixed bonds)
    #4) 

    print("\n--- Create Forcefield ---")
    # --- collect all forces check that all have same cutoff ---
    PairGaussianPotentials = [p for p in Sys.ForceField if p.Names[0]=='ljg']
    ExternalPotentials = [p for p in Sys.ForceField if p.Names[0]=='external_sinusoid']
    BondedPotentials = [p for p in Sys.ForceField if p.Names[0] in ['angle','anglespline','bond','torsion','torsionspline'] ]
    HasEwald=False
    HasEwaldSmearCorrection=False
    SmearedEwaldPotentials = []
    for P in Sys.ForceField:
        if isinstance(P, sim.potential.Ewald) or isinstance(P, sim.potential.ewald.Ewald):
            HasEwald = True
            EwaldPotential = P
        if isinstance(P, sim.potential.smearedcoulombEwaldCorrection.SmearedCoulombEwCorr):
            HasEwaldSmearCorrection = True
            SmearedEwaldPotentials.append(P)
    print('HasEwald {}'.format(HasEwald))
    print('HasEwaldSmearCorrection {}'.format(HasEwaldSmearCorrection))
    #=========================
    #create HarmonicBondForce#
    #=========================
    print("---> Harmonic bond")
    #--- make sure that all bonded site types have a harmonic potential or a rigid bond
    STypes = Sys.World.SiteTypes
    NSTypes = len(STypes)
    BondMatrix = np.zeros((NSTypes, NSTypes), dtype = bool)
    #--- first mark the rigid bonds
    BondOrdMatrix = Sys.World.GetBondOrdMatrix(ShowRigid=True)
    BondMatrix[BondOrdMatrix == -2] = True
    print("Site-Type Bond Matrix: {0}x{0}".format(len(BondMatrix)))
    print(BondMatrix)
    #--- now find the ones with a potential
    for P in Sys.ForceField:
        if P.Type == ptypes.PairPotential and P.Filter.Bonded:
            for (SType1, SType2) in P.Filter.Select(STypes):
                BondMatrix[SType1.SID, SType2.SID] = True
                BondMatrix[SType2.SID, SType1.SID] = True 
    print(np.array(BondMatrix).astype(int))
    #--- now run through all bonds and ensure each had a potential
    for (SType1, SType2) in sim.atomselect.BondPairs.Select(STypes):
        ID1, ID2 = SType1.SID, SType2.SID 
        if not BondMatrix[ID1, ID2]:
            print("... bond matrix test failed for sites {},{}".format(SType1.Name,SType2.Name))
            raise OMMError("Could not find a harmonic or rigid bond between site types {:d} and {:d}".format(ID1, ID2))


    #--- now iterate through all (pair) bonds and add appropriate potential
    #check that only one possible potential applies
    #check that harmonic
    bondedForce = openmm.HarmonicBondForce()

    for ib,b in enumerate(top.bonds()):
        #print(b)
        #Find the filters that work on the bond. Should have exactly one
        ai = b[0].index
        aj = b[1].index
        atomsInThisBond = [Sys.Atom[ai],Sys.Atom[aj]]
        applicablePotentials = [p for p in BondedPotentials if len(list( p.Filter.Select(atomsInThisBond) ))==1]
        #applicablePotentials = []
        #for p in BondedPotentials:
        #    if len( list(p.Filter.Select([Sys.Atom[ai],Sys.Atom[aj]]) ) ) == 1: #should have exactly one bond match
        #        applicablePotentials.append(p)
        #        print(p)

        if len(applicablePotentials) == 0:
            if constraintlengths[ib] == 0:
                raise OMMError("no bonded potential nor constraint found for bond {}".format(b))
        elif len(applicablePotentials) > 1:
            raise OMMError("more than one bonded potential defined for {}, ambiguous".format(b))
       
        #Add the bonded potential
        if len(applicablePotentials) > 0:
            thisP = applicablePotentials[0]
            if thisP.Names[0] == 'bond':
                bondedForce.addBond(ai,aj, thisP.Dist0[0]*sigma, 2.0*thisP.FConst[0]*epsilon/sigma/sigma)      
            else:
                raise OMMError("only harmonic bond implemented right now")

    system.addForce(bondedForce)

    #===================================================
    #create custom nonbonded force: Gaussian + ShortRange (Including ewald smearing correction)+ ExtField#
    #===================================================
    #--- iterate over all atom pair types ---
    #for aname1 in atomNameMap:
    #    for aname2 in atomNameMap:
    #        AType1 = simAtomTypeMap[aname1]
    #        AType2 = simAtomTypeMap[aname2]
    #        #AType1 = [atype for atype in Sys.World.AtomTypes if atype.Name==aname1]
    #        #AType2 = [atype for atype in Sys.World.AtomTypes if atype.Name==aname2]
    #        print([AType1,AType2])
    print("---> Nonbonded Force")
    #--- find global cutoff ---
    GlobalCut = 0
    if np.any(Sys.BoxL > 0):
        MaxCut = 0.5 * np.min(Sys.BoxL[Sys.BoxL > 0])
    else:
        MaxCut = 1.e300
    for P in Sys.ForceField:
        if not P.Type == ptypes.PairPotential: continue
        if not P.Cut is None and P.Cut > MaxCut:
            raise OMMError("Cutoff for potential %s greater (%11.4e) than half min box length (%11.4e)." % (P.Name, P.Cut, MaxCut))
        if GlobalCut < P.Cut:
            GlobalCut = P.Cut
    print("Global cutoff (largest specified cutoff): {}".format(GlobalCut))


    #--- first find nonbonded interactions and smeared Coulomb Ewald corrections, and tabulate accumulate their parameters into tables ---
    HasPair = False
    NTypes = len(Sys.World.AtomTypes) #Assume that sim-Sys defined wiht filters using AtomTypes, not Site Types
    PIndex = [] #2D array indexed by atom types; for each atom pair give list of pair potentials (and an index pointing to the pair) involving said atom pair
    PSmEwCorrIndex = [] #2D array of ewald smearing corrections
    for ID1 in range(NTypes):
        PIndex.append([])
        PSmEwCorrIndex.append([])
        for ID2 in range(NTypes):#range(0, ID1+1):
            PIndex[-1].append([])
            PSmEwCorrIndex[-1].append([])

    for (i, P) in enumerate(Sys.ForceField):
        if isinstance(P, sim.potential.Ewald): continue
                 
        if P.Type == ptypes.PairPotential and not P.Filter.Bonded:
            #select atom types that work for this potential -- why all of a sudden becomes site types? i.e. why can't use Select on AtomTypes?
            HasPair = True
            for (SType1, SType2) in P.Filter.Select(Sys.World.SiteTypes):
                ID1, ID2 = SType1.AID, SType2.AID
                if ID1 < ID2:
                    ID1, ID2 = ID2, ID1
                TypeInd = P.GetTypeInd(SType1, SType2) #for charges, since potential can apply differently (i.e. diff. q1q2) to diff. species pairs
                thisP = (P, TypeInd)
                #add potential to running lists for right bond types
                if isinstance(P, sim.potential.smearedcoulombEwaldCorrection.SmearedCoulombEwCorr): #smeared coulomb correction can't be accumulated for accuracy concerns
                    workingIndex = PSmEwCorrIndex
                else:
                    workingIndex = PIndex

                if not thisP in workingIndex[ID1][ID2]:
                    workingIndex[ID1][ID2].append(thisP)
                    if ID1!=ID2: workingIndex[ID2][ID1].append(thisP)
    

    if UseTabulated:
        knots = np.zeros([NTypes,NTypes,int(NPairPotentialKnots)])
        SPC0 = np.zeros([NTypes,NTypes,NPairPotentialKnots])
        SPC1 = np.zeros([NTypes,NTypes,NPairPotentialKnots])
        SPC2 = np.zeros([NTypes,NTypes,NPairPotentialKnots])
        SPC3 = np.zeros([NTypes,NTypes,NPairPotentialKnots])
            
        print("Using tabulation to aggregate interaction between every pair of atom types")
        for ID1 in range(NTypes):
            for ID2 in range(NTypes):
                PList = PIndex[ID1][ID2]
                rs,us,cut = TabulatePair(Sys, PList, GlobalCut) #get function values to fit spline to
                if ID1==0 and ID2==0:
                    tabulatedpotentials = np.zeros([NTypes,NTypes,int(NPairPotentialBins+1)])
                tabulatedpotentials[ID1,ID2,:] = us
                ks,spcs = GetKnots(GlobalCut,NPairPotentialKnots,rs,us)
                knots[ID1,ID2,:] = ks
                SPC0[ID1,ID2,:] = spcs[0,:]
                SPC1[ID1,ID2,:] = spcs[1,:]
                SPC2[ID1,ID2,:] = spcs[2,:]
                SPC3[ID1,ID2,:] = spcs[3,:]
    else:
        print("!!! Assuming only a single LJGaussian per pair")
        # later, can accumulate the potentials to allow for arbitrary interactions.
        LJGIndex = []
        for ID1 in range(NTypes):
            LJGIndex.append([])
            for ID2 in range(NTypes):#range(0, ID1+1):
                LJGIndex[-1].append([])

        LJGcutoffs = []
        nspec = len(Sys.World.AtomTypes)
        epsmatrix = np.zeros([nspec,nspec])
        sigmatrix = np.ones([nspec,nspec])
        Bmatrix = np.zeros([nspec,nspec])
        kappamatrix = np.zeros([nspec,nspec])
        dist0matrix = np.zeros([nspec,nspec])
        for ID1 in range(NTypes):
            for ID2 in range(NTypes):
                #if len(PIndex[ID1][ID2])!=1:
                #    print("Error: Unsupported number of interactions for atom pair ({},{})".format(ID1,ID2))
                LJGIndex[ID1][ID2] = [ p[0] for p in PIndex[ID1][ID2] if p[0].Names[0] == 'ljg' ]
                if len(LJGIndex[ID1][ID2]) == 0:
                    print("atom type pair ({},{}): No LJ interaction defined, continuing".format(ID1,ID2))
                    continue
                if len(LJGIndex[ID1][ID2]) > 1:
                    print("Error: Unsupported number of LJ-Gaussian interactions for atom pair ({},{})".format(ID1,ID2))
                if len(LJGIndex[ID1][ID2]) != len(PIndex[ID1][ID2]):
                    print("Error: Extra potentials other than LJ-Gaussian detected for atom pair ({},{}); IGNORING THEM".format(ID1,ID2))

                #Store the interaction Matrix
                LJGcutoffs.extend( [p.Cut for p in LJGIndex[ID1][ID2]] )
                print("atom type pair ({},{}): {}".format(ID1,ID2,LJGIndex[ID1][ID2]))
                epsmatrix[ID1,ID2] = LJGIndex[ID1][ID2][0].Epsilon[0]*epsilon.value_in_unit(unit.kilojoule/unit.mole)
                #epsmatrix[ID1,ID2] = LJGIndex[ID1][ID2][0].Epsilon[0]
                sigmatrix[ID1,ID2] = LJGIndex[ID1][ID2][0].Sigma[0]
                Bmatrix[ID1,ID2]   = LJGIndex[ID1][ID2][0].B[0]*epsilon.value_in_unit(unit.kilojoule/unit.mole)
                #Bmatrix[ID1,ID2]   = LJGIndex[ID1][ID2][0].B[0]
                kappamatrix[ID1,ID2] = LJGIndex[ID1][ID2][0].Kappa[0]
                dist0matrix[ID1,ID2] = LJGIndex[ID1][ID2][0].Dist0[0]
        print("Epsilon")
        print(epsmatrix)
        print("Sigma")
        print(sigmatrix)
        print("B")
        print(Bmatrix)
        print("Kappa")
        print(kappamatrix)
        print("Dist0")
        print(dist0matrix)
    #End ifelse for nonbonded pair potential type
    
    #Get Smeared Coulomb Matrix
    if HasEwaldSmearCorrection:
        print("!!! If a pair of species has more than one smeared coulomb correction, using only the last one defined!")
        # later, can accumulate the potentials to allow for arbitrary interactions.
        SmEwCorrIndex = []
        for ID1 in range(NTypes):
            SmEwCorrIndex.append([])
            for ID2 in range(NTypes):#range(0, ID1+1):
                SmEwCorrIndex[-1].append([])

        SmEwCorrCutoffs = []
        nspec = len(Sys.World.AtomTypes)
        aBornMatrix = np.ones([nspec,nspec])
        smPrefactorMatrix = np.zeros([nspec,nspec])
        for ID1 in range(NTypes):
            for ID2 in range(NTypes):
                #if len(PIndex[ID1][ID2])!=1:
                #    print("Error: Unsupported number of interactions for atom pair ({},{})".format(ID1,ID2))
                SmEwCorrIndex[ID1][ID2] = [ p[0] for p in PSmEwCorrIndex[ID1][ID2] if p[0].Names[0] == 'smearedcoulombEwCorr' ]
                if len(SmEwCorrIndex[ID1][ID2]) == 0:
                    print("atom type pair ({},{}): No smeared Ewald Correction interaction defined, continuing".format(ID1,ID2))
                    continue
                if len(SmEwCorrIndex[ID1][ID2]) > 1:
                    print("Error: Unsupported number of smeared Ewald Correction interactions for atom pair ({},{}), only using last one".format(ID1,ID2))
                if len(SmEwCorrIndex[ID1][ID2]) != len(PIndex[ID1][ID2]):
                    print("Error: Extra potentials other than smeared Ewald Correction detected for atom pair ({},{}) were accidentally tabulated together, IGNORING THEM".format(ID1,ID2))

                #Store the interaction Matrix
                SmEwCorrCutoffs.extend( [p.Cut for p in SmEwCorrIndex[ID1][ID2]] )
                print("atom type pair ({},{}): {}".format(ID1,ID2,SmEwCorrIndex[ID1][ID2]))
                smPrefactorMatrix[ID1,ID2] = SmEwCorrIndex[ID1][ID2][-1].Coef[0]
                aBornMatrix[ID1,ID2] = SmEwCorrIndex[ID1][ID2][-1].BornA[0]
        print("Born Radii")
        print(aBornMatrix)
        print("Smearing Interaction Prefactor Matrix: check that same as Ewald")
        print(smPrefactorMatrix)
    #end smeared coulomb matrix parameter accumulation


    #--- iterate over atom types to obtain their external potential fields ---#
    #currently only allow one external field to apply to each atom type. Later on can use tabulated to sum up arbitrary potentials
    PExtIndex = []
    for ID1 in range(NTypes):
        PExtIndex.append([])

    for (ip,P) in enumerate(ExternalPotentials):
        if P.Type == ptypes.FieldPotential:
            for SType1 in P.Filter.Select(Sys.World.SiteTypes):
                if type(SType1) is list: #for some reason, if only one site type, Filter.Select returns a list instead of just that one item!
                    SType1 = SType1[0]
                ID1 = SType1.AID
                if P not in PExtIndex[ID1]:
                    if P.Names[0] == 'external_sinusoid':
                        PExtIndex[ID1].append(P)
                    else:
                        print("Error: unsupported external potential {}".format(P.Names[0]))

    #check that only one external potential per species for now
    #later, can accumulate external potential in a tabulated form
    #can also consider using InteractionGroups(), and adding a diff potential for every filter pair set, just like LAMMPS or sim! Then need to program in custom functional forms to add.
    #  one possible problem is that InteractionGroups() don't use neighbor lists! https://github.com/pandegroup/openmm/issues/1765
    for ID1 in range(NTypes):
        tmp = [p for p in PExtIndex[ID1] if p.Names[0] == 'external_sinusoid']
        if len(tmp) > 1:
            print("Error: Unsupported # of external sinusoids acting on atom {}".format(ID1))
            print(tmp)


    #--- create the nonbonded potential ---
    nonbondedMethod = 2
    if UseTabulated:
        if HasPair: 
            print("Nonbonded interactions detected, adding using tabulated potential")
            
            # == define spline function ==
            nknots = np.ones(NTypes*NTypes) * NPairPotentialKnots
           
            energy_function = 'SPC0(type1,type2,k) + (x-k)*( SPC1(type1,type2,k) + (x-k)*( SPC2(type1,type2,k) + (x-k)*SPC3(type1,type2,k) ) );'
            energy_function += 'k = min(floor(x), nknot(type1,type2)-1);'
            energy_function += 'x=r*nknot(type1,type2)/rcut;'
            fspline = openmm.CustomNonbondedForce(energy_function);
            fspline.addGlobalParameter('rcut',GlobalCut)
            fspline.addPerParticleParameter('type')
            #fspline.setCutoffDistance( fnb.getCutoffDistance() )
            #fspline.setNonbondedMethod(min(fnb.getNonbondedMethod(), 2))
            fspline.setCutoffDistance(GlobalCut) 
            fspline.setNonbondedMethod(nonbondedMethod)

            fspline.addTabulatedFunction('SPC0', openmm.Discrete3DFunction(NTypes,NTypes,NPairPotentialKnots, SPC0.ravel(order='F')))
            fspline.addTabulatedFunction('SPC1', openmm.Discrete3DFunction(NTypes,NTypes,NPairPotentialKnots, SPC1.ravel(order='F')))
            fspline.addTabulatedFunction('SPC2', openmm.Discrete3DFunction(NTypes,NTypes,NPairPotentialKnots, SPC2.ravel(order='F')))
            fspline.addTabulatedFunction('SPC3', openmm.Discrete3DFunction(NTypes,NTypes,NPairPotentialKnots, SPC3.ravel(order='F')))
            fspline.addTabulatedFunction('nknot', openmm.Discrete2DFunction(NTypes,NTypes, nknots.ravel(order='F')))

            # == add types ==
            for atom in top.atoms():
                fspline.addParticle( [atomTypeIndex[atom.name]] )
            
            system.addForce(fspline)
            

            """#Using OpenMM's own tabulation -- slow if want more points (for accuracy), and usually only accurate to ~0.001
            energy_function = "Tabulated(type1,type2,r);"
            fcnb = openmm.CustomNonbondedForce(energy_function)
            fcnb.addPerParticleParameter('type')
            fcnb.setCutoffDistance( GlobalCut )
            fcnb.setNonbondedMethod( nonbondedMethod ) 

            fcnb.addTabulatedFunction('Tabulated',openmm.Continuous3DFunction(NTypes,NTypes,int(NPairPotentialBins+1),tabulatedpotentials.ravel(order='F'), 0,NTypes-1,0,NTypes-1,0,GlobalCut ))
            for atom in top.atoms():
                fcnb.addParticle( [atomTypeIndex[atom.name]] )
            system.addForce(fcnb)
            """ 
    else:
        if len(LJGcutoffs) > 0:
            print("Nonbonded interactions detected, adding. Choosing to assume only one Gaussian interaction per pair")
            if Verbose: print("Detected cutoffs: {}".format(LJGcutoffs))
            if all(cut == LJGcutoffs[0] for cut in LJGcutoffs):
                GlobalCut    = LJGcutoffs[0]
                print("using global cutoff: {}".format(GlobalCut))
            else:
                GlobalCut    = min(LJGcutoffs)
                print("Error: not all cutoffs the same! Using the smallest value: {}".format(GlobalCut))

            energy_function =  'LJ + Gaussian - CutoffShift;'
            energy_function += 'LJ = 0;'
            #energy_function += 'LJ = 4*eps(type1,type2)*(rinv12 - rinv6);'
            energy_function += 'Gaussian = B(type1,type2)*exp(-kappa(type1,type2)*(r - dist0(type1,type2))^2);'
            #energy_function += 'rinv12 = rinv6^2;'
            #energy_function += 'rinv6 = (sig(type1,type2)/r)^6;'
            energy_function += 'CutoffShift = 4*eps(type1,type2)*( (sig(type1,type2)/{0})^12 - (sig(type1,type2)/{0})^6 ) + B(type1,type2)*exp(-kappa(type1,type2)*({0}-dist0(type1,type2))^2);'.format(GlobalCut)

            fcnb = openmm.CustomNonbondedForce(energy_function)
            fcnb.addPerParticleParameter('type')
            fcnb.setCutoffDistance( GlobalCut )
            fcnb.setNonbondedMethod( nonbondedMethod ) #2 is cutoff non periodic

            fcnb.addTabulatedFunction('eps', openmm.Discrete2DFunction(nspec,nspec,epsmatrix.ravel(order='F')) )
            fcnb.addTabulatedFunction('sig', openmm.Discrete2DFunction(nspec,nspec,sigmatrix.ravel(order='F')) )
            fcnb.addTabulatedFunction('B', openmm.Discrete2DFunction(nspec,nspec,Bmatrix.ravel(order='F')) )
            fcnb.addTabulatedFunction('kappa', openmm.Discrete2DFunction(nspec,nspec,kappamatrix.ravel(order='F')) )
            fcnb.addTabulatedFunction('dist0', openmm.Discrete2DFunction(nspec,nspec,dist0matrix.ravel(order='F')) )
            for atom in top.atoms():
                fcnb.addParticle( [atomTypeIndex[atom.name]] )
            system.addForce(fcnb)

        #end ifUseTabulated

   

    #--- create the external potential ---
    print("---> External Force")
    direction=['x','y','z']

    fExts=[]
    for ID1 in range(NTypes):
        if PExtIndex[ID1]:
            thisP = PExtIndex[ID1][0]
            external={"planeLoc":thisP.PlaneLoc, "ax":thisP.PlaneAxis, "U":thisP.UConst[0]*epsilon.value_in_unit(unit.kilojoule/unit.mole), "NPeriod":thisP.NPeriods[0]}

            print("Adding external potential for atom type {} in the {} direction".format(ID1,direction[external["ax"]]))
            energy_function = 'U*sin(2*{pi}*NPeriod*(r-{r0})/{L}); r={axis};'.format(pi=np.pi, L=box_edge[external["ax"]], r0=external["planeLoc"], axis=direction[external["ax"]])
            fExts.append( openmm.CustomExternalForce(energy_function) )
            fExts[-1].addGlobalParameter("U", external["U"])
            fExts[-1].addGlobalParameter("NPeriod", external["NPeriod"])

            for ia,atom in enumerate(top.atoms()):
                if atomTypeIndex[atom.name] == ID1:
                    #print('adding {} to external force'.format(ID1))
                    fExts[-1].addParticle( ia,[] )
            
            system.addForce(fExts[-1])
         
        else:
            print("No external force detected for atom type {}, continuing".format(ID1))
            continue

    for f in fExts:
        print("External potential with amplitude U={}, NPeriod={}".format(f.getGlobalParameterDefaultValue(0), f.getGlobalParameterDefaultValue(1)))

    #=================================
    #Create Nonbonded force for Ewald#
    #=================================
    if HasEwald:
        print('Adding Ewald')
        nbfmethod = openmm.NonbondedForce.PME 
        ewald_tolerance = 1e-5
        print("To implement in OMM, unit charge is now {}".format(qFactor))
        chargeScale = qFactor * EwaldPotential.Coef[0]**0.5

        nbf = openmm.NonbondedForce()
        nbf.setCutoffDistance( GlobalCut )
        nbf.setEwaldErrorTolerance( ewald_tolerance )
        nbf.setNonbondedMethod( nbfmethod )
        nbf.setUseDispersionCorrection(False)
        nbf.setUseSwitchingFunction(False)

        for (i,Atom) in enumerate(Sys.Atom):
            charge = Atom.Charge * chargeScale #In dimensionless, EwaldPotential.Coef is typically 1, and usually change relative strength via temperature. But can also scale the coef, which then acts as lB in the unit length
            LJsigma = 1.0
            LJepsilon = 0.0
            nbf.addParticle(charge, LJsigma, LJepsilon)
        
        system.addForce(nbf)
    #End HasEwald
    
    if HasEwaldSmearCorrection and HasEwald: #Need explicit analytical portion to cancel out 1/r; tabulation can't capture such a steep potential
        #Coef = SmearedEwaldPotentials[0].Coef[0]
        Coef = EwaldPotential.Coef[0]
        if np.all( np.logical_or( smPrefactorMatrix == 0., smPrefactorMatrix == Coef ) ):
            print("Smearing potential coefficients correctly set up")
        else:
            print("WARNING! smearing potential coefficients suspicious, do not match ewald coefficient, using the ewald coefficient.")
        print( "Ewald coefficient is: {}".format(Coef))
        print( "Matrix of smearing correction coefficients is:" )
        print(smPrefactorMatrix)

        energy_function = 'coef*q1*q2 * ( (erf(factor*r) - 1)/r - shift );'
        energy_function += 'shift = (erf(factor*rcut) -1)/rcut;'
        energy_function += 'factor = sqrt({:f})/2/aborn(type1,type2);'.format(np.pi)
        energy_function += 'coef = {:f};'.format(Coef)
        energy_function += 'rcut = {:f};'.format(GlobalCut)
        fcnb = openmm.CustomNonbondedForce(energy_function)
        fcnb.addPerParticleParameter('type')
        fcnb.addPerParticleParameter('q')
        fcnb.setCutoffDistance( GlobalCut )
        fcnb.setNonbondedMethod( openmm.NonbondedForce.CutoffPeriodic )
        fcnb.addTabulatedFunction('aborn',openmm.Discrete2DFunction(nspec, nspec, aBornMatrix.ravel(order='F')))
        for ia,atom in enumerate(top.atoms()):
            q = Sys.Atom[ia].Charge
#            print("atom {} has charge {}".format(ia,q))
            fcnb.addParticle( [atomTypeIndex[atom.name], q] )
        system.addForce(fcnb)

    return  system,top,mdtrajtop


def CreateOMMSimulation(Sys,system,top,Prefix="",chkfile='checkpnt.chk',Verbose=False):
    # --- integration options ---
    # TODO: NPT
    reduced_timestep = Sys.Int.Method.TimeStep#0.01
    reduced_temp = Sys.TempSet#1
    reduced_Tdamp = Sys.Int.Method.LangevinGamma#1 #time units
    dt = reduced_timestep * tau
    temperature = reduced_temp * epsilon/kB/N_av
    friction = 1/(reduced_Tdamp) /tau
    
    pressure = Sys.PresSet
    if Sys.PresSet > 0:
        useNPT = True
        useNVT = False
    else:
        useNPT = False
        useNVT = True
    """
    reduced_pressure = 1
    reduced_Pdamp = 0.1 #time units
    pressure = reduced_pressure * epsilon/(sigma**3) / N_av
    barostatInterval = int(reduced_Pdamp/reduced_timestep)
    """

    #===========================
    ## Prepare the Simulation ##
    #===========================
    print("=== Preparing Simulation ===")
    if useNPT:
        pressure = pressure * epsilon/N_av/sigma/sigma/sigma #convert from unitless to OMM units
        barostatInterval = 25 #in units of time steps. 25 is OpenMM default
        system.addForce(openmm.MonteCarloBarostat(pressure, temperature, barostatInterval))
        print("Added MC Barostat with P {} (eps/sig^3), T {}, freq {}".format(
            6.02214179e23 * pressure.value_in_unit(unit.kilojoules/unit.nanometer**3),temperature,barostatInterval))
        print("In OMM units, is P={}'bar'".format(pressure.value_in_unit(unit.bar)))

    integrator = openmm.LangevinIntegrator(temperature, friction, dt)
    if system.getNumConstraints() > 0:
        print("Applying bond constraints before starting")
        integrator.setConstraintTolerance(constraintTolerance)

    # TODO: simulation platform
    #platform = openmm.Platform.getPlatformByName('CUDA')
    #platformProperties = {'Precision': 'mixed'}
    #platform = openmm.Platform.getPlatformByName('CPU')
    #platform = openmm.Platform.getPlatformByName('Reference')
    #simulation = app.Simulation(top,system, integrator, platform)
    #simulation = app.Simulation(top,system, integrator, platform,platformProperties)
    if device == -1:
        print("Automatically choosing platform")
        simulation = app.Simulation(top,system,integrator)
    else:
        print("Manually setting platform: {} and device: {}".format(platformName,device))
        platform = openmm.Platform.getPlatformByName(platformName)
        if platformName == 'CUDA':
            platform.setPropertyDefaultValue("CudaDeviceIndex", str(device))
        elif platformName == 'OpenCL':
            platform.setPropertyDefaultValue("OpenCLDeviceIndex", str(device))
        else:
            platform = openmm.Platform.getPlatformByName('CPU')
        simulation = app.Simulation(top,system, integrator, platform)

    chkfile = Prefix + chkfile
    simulation.reporters.append(app.checkpointreporter.CheckpointReporter(chkfile,10000))


    # --- done ---
    simOptions = {'dt':dt, 'temp':temperature, 'fric':friction}
    print(simOptions)
    return simulation,simOptions 


def DCDtoLAMMPS(TrajFile,topfile):
    '''Converts DCD to LAMMPS. Right now need to divide by 10 because of default mdtraj behavior. Should modify to set mdtraj's default unit to dimensionless.
    '''
    import mdtraj
    trj = mdtraj.load(TrajFile, top=topfile)
    trj.xyz/=10.0
    trj.unitcell_lengths/=10.0
    outname = '.'.join( TrajFile.split('.')[:-1] ) + '.lammpstrj'
    trj.save(outname)

    return outname


def RunOpenMM(Sys, Prefix = "", Verbose = False, NStepsMin = 0, NStepsEquil = 10000, NStepsProd = 10000, WriteFreq = 10000, TrajFile = "trj.dcd"):
    "Run OpenMM. Returns OMMFiles (initial pdb, equilibration pdb, final pdb, dcd, thermo log, output log, simulation preamble, returncode)"
    OpenMMFiles = []
    # --- Get OMM System and Simulation ---
    system,top,mdtrajtop = CreateOMMSys(Sys,Verbose)
    simulation,simOptions = CreateOMMSimulation(Sys,system,top,Prefix,Verbose=Verbose)

    # --- simulation options ---
    TrajFile = Prefix + TrajFile
    EqLogFile = Prefix + "eqlog.txt"
    LogFile = Prefix + "prodlog.txt"
    OpenMMFiles.append(LogFile)

    # --- Init ---
    print("=== Initializing Simulation ===")
    Pos = Sys.Pos.copy()
    nImage = sim.geom.MinimageN(Pos, Sys.BoxL)
    Pos = sim.geom.Minimage(Pos, Sys.BoxL)
    
    #need to make molecules whole
    Lx,Ly,Lz = Sys.BoxL[0], Sys.BoxL[1], Sys.BoxL[2]
    nframes = Pos.shape[0]
    unitcell_lengths = np.array([Lx,Ly,Lz])
    print('unitcell_lengths {}'.format(unitcell_lengths))
    unitcell_angles = np.array([90., 90., 90.])
    mdt_traj = mdtraj.Trajectory(np.array(Pos), topology = mdtrajtop, unitcell_lengths = unitcell_lengths, unitcell_angles = unitcell_angles)
    bondlist = np.array([ [b[0].index,b[1].index] for b in top.bonds() ], dtype=np.int32) #should be sorted, as long as the bonds come out sorted from the sim Sys object
    if len(bondlist) > 0:
        print(">0 bonds, making trajectory whole to be safe")
        wholetraj = mdt_traj.make_molecules_whole(inplace=False,sorted_bonds=bondlist)
        simulation.context.setPositions(wholetraj.xyz[0])#Pos)
    else:
        simulation.context.setPositions(Pos)

    initstatefile = Prefix + 'output.xml'
    simulation.saveState(initstatefile)
    OpenMMFiles.append(initstatefile)
    initialpdb = Prefix + "initial.pdb"
    initial_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    app.PDBFile.writeModel(simulation.topology, initial_positions, open(initialpdb,'w'))
    OpenMMFiles.append(initialpdb)

    #TODO: apply constraint tolerance
    if system.getNumConstraints() > 0:
    #if True:
        simulation.context.applyConstraints(constraintTolerance)
        constrainedpdb = Prefix + "constrained.pdb"
        constrained_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
        app.PDBFile.writeModel(simulation.topology, constrained_positions, open(constrainedpdb,'w'))
        OpenMMFiles.append(constrainedpdb)


    #temporary test:
    #for x in [0.01, 0.5, 1., 2., 3., 4.]:
    #    simulation.context.setPositions( [[0.,0.,0.],[x,0.,0.]] )
    #    state = simulation.context.getState(getEnergy=True)
    #    print( "r: {}, energy: {}".format(x, state.getPotentialEnergy()) )

    if NStepsMin > 0:
    #if False:
        print("=== Running energy minimization ===")
        minimizefile = Prefix + "minimized.pdb"
        simulation.minimizeEnergy(maxIterations=NStepsMin)
        simulation.context.setVelocitiesToTemperature(simOptions["temp"]*3)
        minimized_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
        app.PDBFile.writeModel(simulation.topology, minimized_positions, open(minimizefile,'w'))
        OpenMMFiles.append(minimizefile)

    # --- Equilibrate ---
    print('=== Equilibrating ===')
    simulation.reporters.append(app.StateDataReporter(sys.stdout, WriteFreq*100, step=True, potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
    simulation.reporters.append(app.StateDataReporter(EqLogFile, WriteFreq, step=True, potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
    print("Progress will be reported every {} steps".format(WriteFreq))
    simulation.step(int(NStepsEquil))
    equilibratefile = Prefix + "equilibrated.pdb"
    equilibrated_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    app.PDBFile.writeModel(simulation.topology, equilibrated_positions, open(equilibratefile,'w'))
    OpenMMFiles.append(equilibratefile)
    

    # --- Production ---
    print('=== Production {} steps total, {} tau total ({}) ==='.format(NStepsProd, NStepsProd*Sys.Int.Method.TimeStep, NStepsProd*Sys.Int.Method.TimeStep*tau ))
    simulation.reporters.pop()
    simulation.reporters.append(app.StateDataReporter(LogFile, WriteFreq, step=True, potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True, density=True, speed=True, separator='\t'))
    simulation.reporters.append(mdtraj.reporters.DCDReporter(TrajFile,WriteFreq))
    #simulation.reporters.append(app.checkpointreporter.CheckpointReporter('checkpnt.chk', 5000))
    start = time.time()
    simulation.step(int(NStepsProd))
#    simulation.step(NStepsProd)
    end = time.time()
    print("=== Finished production in {} seconds".format(end-start))

    # --- Returns ---
    #TODO: validate output for success
    returncode = 0
    return OpenMMFiles, TrajFile, returncode


def MakeOpenMMTraj(Sys, DelTempFiles = False, Prefix = "", Verbose = False, *args, **kwargs):
    """Runs OpenMM and converts into LAMMPS trajectory object that sim can read."""    
    import sim.traj
    import os
    
    OpenMMFiles, DCDFile, returncode = RunOpenMM(Sys, *args, Prefix = Prefix, Verbose = Verbose, **kwargs)

    #convert DCD to LAMMPS
    initialpdb = [filename for filename in OpenMMFiles if "initial" in filename][0]
    LogFile = [filename for filename in OpenMMFiles if "prodlog" in filename][0]
    if Verbose:
        print("=== Temporary files from openMM simulation: ===")
        print(OpenMMFiles)
        print(initialpdb)
        print(LogFile)

    TrajFile = DCDtoLAMMPS(DCDFile, topfile=initialpdb)

    #usual stuff --> need to create new trajectory function, to parse out thermo info
    #Traj = sim.traj.Lammps(TrajFile,LogFile)
    
    Traj = sim.traj.Lammps(TrajFile,LogFile,OMMLog=True)
    
    if DelTempFiles:
        for fn in OpenMMFiles:
            if os.path.isfile(fn):
                print('deleting {}'.format(fn))
                os.remove(fn)

    sys.stdout.flush()
    return Traj, TrajFile, DCDFile

def MakeRerunOpenMMTraj(Sys, ElecTrajFile,ElecTrajFileDCD, DelTempFiles = False, Prefix = "", Verbose = False, *args, **kwargs):
    """Reprocess trajectory using forcefield of Sys"""
    import sim.traj
    import os

    LogFile = RerunOpenMM(Sys, ElecTrajFileDCD, *args, Prefix = Prefix, Verbose = Verbose, **kwargs)

    Traj = sim.traj.Lammps(ElecTrajFile,LogFile,OMMLog=True)

    sys.stdout.flush()
    return Traj

def RerunOpenMM(Sys, ElecTrajFileDCD, Prefix = "", Verbose = False, NStepsMin = 0, NStepsEquil = 10000, NStepsProd = 10000, WriteFreq = 1, TrajFile = "trj.dcd"):
    "Run OpenMM. Returns OMMFiles (initial pdb, equilibration pdb, final pdb, dcd, thermo log, output log, simulation preamble, returncode)"
    OpenMMFiles = []
    # --- Get OMM System and Simulation ---
    system,top,mdtrajtop = CreateOMMSys(Sys,Verbose)
    simulation,simOptions = CreateOMMSimulation(Sys,system,top,Prefix,Verbose=Verbose)

    # --- simulation options ---
    LogFile = Prefix + "prodlog.txt"

    # --- Init ---
    print("=== Initializing Simulation ===")
    Pos = Sys.Pos.copy()
    nImage = sim.geom.MinimageN(Pos, Sys.BoxL)
    Pos = sim.geom.Minimage(Pos, Sys.BoxL)

    #need to make molecules whole
    Lx,Ly,Lz = Sys.BoxL[0], Sys.BoxL[1], Sys.BoxL[2]
    nframes = Pos.shape[0]
    unitcell_lengths = np.array([Lx,Ly,Lz])
    unitcell_angles = np.array([90., 90., 90.])
    mdt_traj = mdtraj.Trajectory(np.array(Pos), topology = mdtrajtop, unitcell_lengths = unitcell_lengths, unitcell_angles = unitcell_angles)
    bondlist = np.array([ [b[0].index,b[1].index] for b in top.bonds() ], dtype=np.int32) #should be sorted, as long as the bonds come out sorted from the sim Sys object
    if len(bondlist) > 0:
        print(">0 bonds, making trajectory whole to be safe")
        wholetraj = mdt_traj.make_molecules_whole(inplace=False,sorted_bonds=bondlist)
        simulation.context.setPositions(wholetraj.xyz[0])#Pos)
    else:
        simulation.context.setPositions(Pos)

    initialpdb = Prefix + "initial.pdb"
    initial_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    app.PDBFile.writeModel(simulation.topology, initial_positions, open(initialpdb,'w'))

    #TODO: apply constraint tolerance
    if system.getNumConstraints() > 0:
        simulation.context.applyConstraints(constraintTolerance)
        constrainedpdb = Prefix + "constrained.pdb"
        constrained_positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
        app.PDBFile.writeModel(simulation.topology, constrained_positions, open(constrainedpdb,'w'))
        OpenMMFiles.append(constrainedpdb)


    # --- Reprocess trajectory ---
    print('=== Reprocess trajectory file {} with {} forcefield'.format(ElecTrajFileDCD, Sys.Name))
    traj = mdtraj.load(ElecTrajFileDCD, top = initialpdb)
    nframes = traj.n_frames
    data = []
    for i in range(nframes):
        context = simulation.context
        Pos = traj.xyz[i]
        context.setPositions(Pos)
        PE = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        Vol = context.getState().getPeriodicBoxVolume().value_in_unit((unit.nanometer)**3)
        data.append([i, PE, Vol])
#    print(data)
    np.savetxt(LogFile, np.array(data), header = """\"Step" "Potential Energy (kJ/mole)" "Box Volume (nm^3)\" """)   
    # --- Returns ---
    return LogFile
