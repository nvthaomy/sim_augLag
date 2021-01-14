import numpy as np
import copy,re,random,os,time



class GenerateInchi(object):
    def __init__(self):
        self.reducednetworklist=None
        self.networklist=None
    
    def GetUniqueInchi(self,reducednetworkID):
        self.s=""
        self.blist=copy.copy(self.reducednetworklist[reducednetworkID])
        self.plist=[]
        self.clist=[[]]
        self.t=min(min(self.blist))
        self.newt=None
        thisinchi=''.join([self.sitelist[x] for x in range( min(min(self.networklist[reducednetworkID]))-1,max(max(self.networklist[reducednetworkID]))-1)  ])## UPDATE the zin fork with this edit
        self.Handle()
        thisinchi+=self.s
        return thisinchi

    def RemoveBond(self):
        if [self.t,self.newt] in self.blist:
            self.blist.remove([self.t,self.newt])
        else:
            self.blist.remove([self.newt,self.t])
        self.t=self.newt
    
    def GetBonds(self):
        thesebonds=[]
        for i in self.blist:
            if self.t in i and i not in self.plist and i[::-1] not in self.plist and not any([i in j for j in self.clist]) and not any([i[::-1] in j for j in self.clist]) :
                thesebonds.append(i)
        return copy.copy(thesebonds)
    
    def Handle(self):
        if self.s is '':self.s="%d"%self.t
        thisblist=self.GetBonds()
        if len(thisblist)==0:
            if len(self.clist[-1])>0:
                thisbond=self.clist[-1].pop(-1)
                self.t=thisbond[-1]
                self.blist.remove(thisbond)
                self.s+=",%d"%self.t
                self.Handle()
            elif any([self.t in i for i in self.plist]):
                thisbond=self.plist.pop(-1)
                self.t=thisbond[-1]
                if thisbond in self.blist:
                    self.blist.remove(thisbond)
                else:
                    self.blist.remove(thisbond[::-1])
                self.s+=")%d"%self.t
                self.clist.pop(-1)
                self.Handle()
            else:
                if len(self.plist)>0:
                    thisbond=self.plist[-1]
                    self.t=thisbond[-1]
                    if thisbond in self.blist:
                        self.blist.remove(thisbond)
                    else:
                        self.blist.remove(thisbond[::-1])
                    self.s+=")%d"%self.t
                    self.plist.pop(-1)
                    self.clist.pop(-1)
                    self.Handle()
                else:
                    return
        if len(thisblist)==1:
            self.newt=list(set(thisblist[0])-set([self.t]))[0]
            self.RemoveBond()
            self.s+="-%d"%self.t
            self.Handle()
        if len(thisblist)>=2:
            self.pbond=[self.t,  list(set(thisblist[-1])-set([self.t]))[0]]
            self.plist.append(self.pbond)
            self.clist.append(thisblist[1:-1])
            self.newt=list(set(thisblist[0])-set([self.t]))[0]
            self.RemoveBond()
            self.s+="(%d"%self.t
            self.Handle()
        return







class InterpretTopology(GenerateInchi):
    def __init__(self):
        GenerateInchi.__init__(self)
        self.molpointers,self.networklist=self.GetDiscreteNetworks(self.bondlist)
        self.reducednetworklist=[]
        for network in self.networklist:
            minid=min(min(network))
            self.reducednetworklist.append( [[i-minid,j-minid] for i,j in network] )
        
        uniquenetworklist=[]
        reducedinchilist=[]
        for ii,xblist in enumerate(self.reducednetworklist):
            thisinchi=self.GetUniqueInchi(ii)
            if thisinchi not in reducedinchilist:
                reducedinchilist.append(thisinchi)
                uniquenetworklist.append(self.networklist[ii])
        
        #make a list of unique monomer types (molecules with no bonds)
        bondparticipants=[]#list starting from 1
        for i,j in self.bondlist:
            if i not in bondparticipants:
                bondparticipants.append(i)
            if j not in bondparticipants:
                bondparticipants.append(j)
        
        monomers=[]
        for ii,i in enumerate(self.sitelist):
            if ii+1 not in bondparticipants and i not in monomers:
                monomers.append(i)
        
        #make a molecule list
        molnumlist=[]
        for ii,site in enumerate(self.sitelist):
            Found=False
            for jj,network in enumerate(self.networklist):
                if any([ii+1 in kk for kk in network]):
                    #record the molecule network number
                    molnumlist.append(('molecule',jj))
                    Found=True
            if not Found:
                for jj,monomer in enumerate(monomers):
                    if site == monomer:
                        molnumlist.append(('monomer',jj))
        
        mollist=[]
        previous=None
        for mtype,num in molnumlist:
            if mtype=='monomer':
                mollist.append((mtype,num))
            elif not (mtype,num) == previous:
                mollist.append((mtype,num))
            previous=(mtype,num)
        
        molidlist=[]
        nuniq=len(reducedinchilist)
        #reduce the mollist down to unique moltypes
        for mtype,num in mollist:
            if mtype=='monomer':
                molidlist.append(nuniq+num)
            else:
                network=self.networklist[num]          
                #get the reduced network identifier
                minid=min(min(network))
                reducednetwork= [[i-minid,j-minid] for i,j in network]
                #get the unique molecule index
                molindex=self.reducednetworklist.index(reducednetwork)
                thisinchi=self.GetUniqueInchi(molindex)
                molidlist.append(reducedinchilist.index(thisinchi))
        
        #respool the uniquenetworklist to make it contiguous
        #this will be used for building the bonding geometry of the sysobject
        respooledbonds=[]
        lastmax=0
        correctorlist=[]
        for network in uniquenetworklist:
            thismin=min(min(network))
            thisrespool=[[i-thismin+lastmax , j-thismin+lastmax] for i,j in network]
            correctorlist.append(lastmax-thismin)
            respooledbonds.append(thisrespool)
            lastmax=max(max(thisrespool))+1

        #build a list of siteids for each network
        #these will be used to reference items in the trajectory frame positions array
        networksiteids=[]
        for network in uniquenetworklist:
            networksiteids.append(np.unique(np.array(network).flatten()).tolist())
        
        #build a list of sitetypes for each network
        #these will be used to build the molecule types for the sys object
        networksitetypes=[]
        for networksites in networksiteids:
            networksitetypes.append([self.sitelist[siteid-1] for siteid in networksites])
        
        #make a list of unique site types
        sitetypes=[]
        for stype in self.sitelist:
            if stype not in sitetypes:
                sitetypes.append(stype)
        
        site_dict={}
        for ii,i in enumerate(sitetypes):
            site_dict[i]=ii
        
        self.site_dict=site_dict
        self.sitetypes=sitetypes
        self.networksitetypes=networksitetypes
        self.networksiteids=networksiteids
        self.molidlist=molidlist
        self.monomers=monomers        
        self.respooledbonds=respooledbonds
        self.uniquenetworklist=uniquenetworklist
        self.correctorlist=correctorlist
    
    def GetDiscreteNetworks(self,bondnetwork):
        cp1=list(bondnetwork)
        uniquenetworks=[]
        while len(cp1)>0:
            bondx=cp1[0]
            thisnetwork=[]
            KeepCrawling=True
            bondlist=[]
            while KeepCrawling:
                FoundSomething=0
                cp2=list(cp1)
                for bondy in cp1:
                    if bondy.count(bondx[0])==1 or bondy.count(bondx[1])==1:
                        cp2.remove(bondy)
                        thisnetwork.append(bondy)
                        bondlist.append(bondy)
                        FoundSomething=1
                cp1=cp2
                if FoundSomething:
                    bondx=bondlist.pop()
                if not FoundSomething:
                    if not len(bondlist)==0:
                        bondx=bondlist.pop()
                    else:
                        KeepCrawling=False
            uniquenetworks.append(thisnetwork)
        pointers=[[min(min(network))] for network in uniquenetworks]
        return pointers, uniquenetworks

    def PrmTopBondInfo(self,filename,QueryWater=0,IgnoreWater=True):
        f=file(filename,'r').read()
        f=f.split('%FLAG')
        prmtop=''
        for block in f:
            if block.__contains__('RESIDUE_LABEL'):
                block=block.split('%FORMAT(20a4)                                                                   \n')[1].rstrip()
                reslist=block.split()
            if block.__contains__('BONDS'):
                block=block.split('%FORMAT(10I8)                                                                   \n')[1]
                prmtop+= block
            if block.__contains__('RESIDUE_POINTER'):
                block=block.split('%FORMAT(10I8)                                                                   \n')[1]
                residue_pointer=[int(i) for i in block.split()]
            if block.__contains__('ATOM_NAME'):
                block=block.split('%FORMAT(20a4)                                                                   \n')[1].rstrip()
                block=block.replace('\n','')
                seq=[]
                for i in range(0,len(block),4):
    	            seq.append(block[i:i+4])
                sequence=[x[0] for x in seq]
        f=prmtop.split()
        f=[int(int(x)/3.+1) for x in f]
        f=[f[i:i+2] for i in range(0,len(f),3)]#this is the bondlist
        TerminalsUnique=True#query_yes_no('Make terminal amino acids unique site types?',default='yes')
        if 'WAT' in reslist and QueryWater:
            IgnoreWater=0# query_yes_no('Ignore water molecules?',default='yes')
        if IgnoreWater and 'WAT' in reslist:	
            a=residue_pointer[:reslist.index('WAT')]
            b=sequence[:residue_pointer[reslist.index('WAT')]-1]
            c=reslist[:reslist.index('WAT')]
            d=[bond for bond in f if not (bond[0]>=residue_pointer[reslist.index('WAT')] or bond[1]>=residue_pointer[reslist.index('WAT')])]
            prmbonds,prmpoint,prmsitelist,peptide = d,a,b,c
        else:
            prmbonds,prmpoint,prmsitelist,peptide = f,residue_pointer,sequence,reslist
        molpointers=self.GetDiscreteNetworks(prmbonds)[0]
        aapoint=[[prmpoint.index(x)] for x in [y[0] for y in molpointers]]
        if TerminalsUnique: 
            #combine any capping groups with terminal amino acids, then uniquely label terminal groups with N or C
            if 'ACE' not in peptide:
                for i in aapoint:
    	            if peptide[i[0]]!='WAT': peptide[i[0]]+='N'
            if 'NME' not in peptide and 'NHE' not in peptide:
                if len(aapoint)>1:
    	            for i in aapoint[1:]:
			            if peptide[i[0]-1]!='WAT': peptide[i[0]-1]+='C'
    	            if IgnoreWater:
    		            peptide[-1]+='C'
                else:
    	            if peptide[i[0]-1]!='WAT': peptide[i[0]-1]+='C'
            while 'ACE' in peptide:
                ind=peptide.index('ACE')
                peptide.pop(ind)
                prmpoint.pop(ind+1)
                peptide[ind]+='N'
            while 'NME' in peptide or 'NHE' in peptide:
                try:
    	            ind=peptide.index('NME')
                except:
    	            ind=peptide.index('NHE')
                peptide.pop(ind)
                prmpoint.pop(ind)
                peptide[ind-1]+='C'
        return prmbonds,prmpoint,prmsitelist,[i[0]+i[1:].lower() for i in peptide],[[i+1,i+1] for i,z in enumerate(prmsitelist)]



##############################################################################################################
   
class CreateCG(InterpretTopology):
    def __init__(self,Prefix='a'):
        self.Prefix=Prefix
        self.element_dict={
        'H':[.25,(.9,.9,.9),1,1],
        'C':[.7,(.2,.2,.2),12,1],
        'O':[.6,(1,0,0),16,1],
        'Cl':[1.0,(0,1,0),35,1],
        'N':[.65,(0,0,1),14,1],
        'S':[1.0,(1,1,0),32,1],
        'Se':[1.15,(1,.5,0),79,1],
        }
        prmtopfilename=self.Prefix+'.prmtop.parm7'#'dar16dimer_375K.prmtop.parm7'
        if not os.path.isfile('.parsed.'+prmtopfilename):
            print 'parsing prmtop.parm7'
            self.bondlist,self.pointers,self.sitelist,self.resseq,self.e=self.PrmTopBondInfo(prmtopfilename,IgnoreWater=False)
            fs=open('.parsed.'+prmtopfilename,'w')
            np.save(fs,(self.bondlist,self.pointers,self.sitelist,self.resseq,self.e))
            fs.close()
        else:
            print 'located parsed prmtop file'
            self.bondlist,self.pointers,self.sitelist,self.resseq,self.e=np.load('.parsed.'+prmtopfilename)
        InterpretTopology.__init__(self)
        
        #build a list of the residue types for each network
        #these will be used when the mapping operation is underway to create a mapping rule for each amino acid type
        self.sequences=[]
        self.seqpointers=[]
        for network in self.uniquenetworklist:
            thisseq=[]
            thispoint=[]
            startid=min(min(network))
            stopid=max(max(network))
            for ii,residuestart in enumerate(self.pointers):
                if residuestart<stopid+1 and residuestart>=startid:
                    thisseq.append(self.resseq[ii])
                    thispoint.append(residuestart)
            self.sequences.append(thisseq)
            self.seqpointers.append(thispoint)
        
        #respool the sequence pointers,and put in start,stop format
        self.respooledseqpointers=[]
        for ii,points in enumerate(self.seqpointers):
            thesepointers=[]
            thislength= len(self.networksitetypes[ii])
            for iii,thispointer in enumerate(points):
                if iii == len(points)-1:
                    thesepointers.append( [thispointer+self.correctorlist[ii], points[0]+self.correctorlist[ii]+thislength-1])
                else:
                    thesepointers.append( [thispointer+self.correctorlist[ii], points[iii+1]+self.correctorlist[ii]-1 ])
            self.respooledseqpointers.append(thesepointers)
        
        self.trajids=[]
        for ii,i in enumerate(self.networksitetypes):
            start=self.seqpointers[ii][0]-1
            length=len(i)
            self.trajids.append( range(start,start+length))
        
        #create a convenient container for the molecule information
        self.molecule=zip(self.networksitetypes,self.respooledbonds,self.sequences,self.respooledseqpointers,self.trajids)
        if not os.path.isfile(self.Prefix+'.mapfile'):
            self.Sys=self.MakeMappingSys()
            import sim.system.visualize as vis
            self.vis=vis
            self.Sys.Vis = self.vis.Visualizer3D(self.Sys,ShowBonds=True,Label=0)
            self.Sys.Vis.Update()
            self.Sys.Vis.Display.autocenter=False
            self.Sys.Vis.Display.center=(0,0,0)
            self.cg_dict={}
            self.cgsite_dict={}
            self.MappingDefinition()
        else:
            print 'located predefined mapfile, making generator'
            namespace={}
            execfile(self.Prefix+'.mapfile',namespace)
            self.cg_dict=namespace['cgdict']
            self.element_dict=namespace['element_dict']
        self.MapfileToGenerator()
        ReadGenerator(self.Prefix).MakeSysScript()
        return
    
    def MappingDefinition(self):
        lastmol=None
        for self.molid,(sitelist,bonds,sequence,sequencepointers,trajids) in enumerate(self.molecule):
            self.hideall()
            for aaid, self.aa in enumerate(sequence):
                if self.aa not in self.cg_dict:
                    self.formula=''
                    atomlabels=[]
                    self.cgspheres=[]
                    self.cgrods=[]
                    cglabels=[]
                    print 'amino acid', self.aa
                    if self.molid!=lastmol:
                        oldtail=None
                    lastmol=self.molid
                    atoms=range(sequencepointers[aaid][0],sequencepointers[aaid][1]+1)
                    self.reducedatoms=[atomj-sequencepointers[0][0] for atomj in atoms]
                    self.showghostmol(self.molid)
                    for atomj in self.reducedatoms:
                        self.showatom(self.molid,atomj)
                    possum =np.zeros(self.Sys.Dim, dtype=float)
                    lastpos=np.zeros(self.Sys.Dim, dtype=float)
                    for atomj in self.reducedatoms:
                        thispos=lastpos+self.Sys.MinImage(self.Sys.Mol[self.molid][atomj].Pos-lastpos)
                        possum+=thispos
                        lastpos=thispos
                    center=possum/len(atoms)
                    self.Sys.Pos-=center
                    for atomj in self.reducedatoms:
                        apos=self.Sys.Mol[self.molid][atomj].Pos
                        minpos=apos-self.Sys.BoxL*np.round(apos/self.Sys.BoxL)
                        self.vscale=np.max(self.Sys.BoxL)
                        atomlabels.append(self.vis.v.label(pos=minpos/self.vscale, text=str(self.reducedatoms.index(atomj)+1)))
                    self.Sys.Vis.Update()
                    self.groups=self.InputGroups()
                    for self.group in self.groups:
                        for ii,atom in enumerate(self.group):
                            atomlabels[atom-1].linecolor=(1,0,0)
                            atomlabels[atom-1].color=(1,0,0)
                        self.Sys.Vis.Update()
                        self.newsitename=self.InputSiteName()
                        self.formula+=self.newsitename
                        for ii,atom in enumerate(self.group):
                            atomlabels[atom-1].linecolor=(1,1,1)
                            atomlabels[atom-1].color=(1,1,1)
                        gcenter,gcolor,grg=self.GetProperties()
                        self.cgspheres.append(self.vis.v.sphere(pos=gcenter/self.vscale,color=gcolor,radius=grg/self.vscale,opacity=0.5))
                    for iii,sphere in enumerate(self.cgspheres):
                        cglabels.append(self.vis.v.label(pos=sphere.pos,text=str(iii+1)))
                    for atomlabel in atomlabels:
                        atomlabel.set_visible(False)
                    self.InputConnectivity()
                    time.sleep(3)
                    for cglabel in cglabels:
                        cglabel.set_visible(False)
                    for sphere in self.cgspheres:
                        sphere.set_visible(False)
                    for rod in self.cgrods:
                        rod.set_visible(False)
                    self.Sys.Vis.Update()
        fs=self.WriteMapFileString()
        with open(self.Prefix+'.mapfile','w') as f:
	        f.write(fs)
    
    ################################################CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    def MapfileToGenerator(self):
        if 'Wat' in self.cg_dict:
            IgnoreWater=False
        else:
            IgnoreWater=True
        self.aapoint=[[self.pointers.index(x)] for x in [y[0] for y in self.molpointers]]
        cgbondlist,cgpointer,cgsitelist,inchimap=self.GetCGInfo()
        ambermap=self.GetAAtoCGMap(self.pointers,inchimap,self.resseq)
        bondsets,anglesets,torsionsets,nonbondsets=self.GetAtomSets(cgsitelist,cgbondlist)
        s=''
        s+='\n\ncgsitelist='+str(cgsitelist)
        s+='\n\ncgbondlist='+str(cgbondlist)
        s+='\n\nambermap='+str(ambermap)
        s+='\n\nbondsets='+str(bondsets)
        s+='\n\nanglesets='+str(anglesets)
        s+='\n\ntorsionsets='+str(torsionsets)
        s+='\n\nnonbondsets='+str(nonbondsets)
        s+="\n\nelement_dict={"
        for k,v in self.element_dict.items():
            s+="""\n"%s" : %s,"""%(str(k),str(v))
        s.rstrip(',')
        s+="\n}\n\n"
        with open(self.Prefix+'_generatorCG.txt','w') as f:
            f.write(s)
    
    def GetAtomSets(self,sitelist,bondlist):
        bondsets=[]
        anglesets=[]
        torsionsets=[]
        nonbondsets=[]
        bonds=self.FindBondPairs(sitelist,bondlist)
        bondsets.extend(['Bond'+i for i in bonds])
        triplesid=self.FindBondTriples(sitelist,bondlist)
        triples=self.RemovePalindromes(set([''.join([sitelist[i-1] for i in trip]) for trip in triplesid]))
        anglesets.extend(['Angle'+i for i in triples])
        torsions=self.FindBondQuartets(sitelist,bondlist,triplesid)
        torsionsets.extend(['Torsion'+i for i in torsions])
        nonbondsets.extend(['NonBond'+i for i in self.FindNonBondPairs(sitelist,bonds,triples,torsions)])
        return (bondsets,anglesets,torsionsets,nonbondsets)    
    
    def FindBondPairs(self,sitelist,bondlist):
        pairs=[]
        for pair in bondlist:
            pairs.append(''.join([sitelist[i-1] for i in pair]))
        return self.RemovePalindromes(set(pairs))
    
    def FindNonBondPairs(self,sitelist,bonds,triples,torsions):
        safelist=[j[0]+j[2] for j in [[i for i in re.findall('[A-Z][^A-Z]*',x)] for x in triples]]
        safelist.extend([j[0]+j[3] for j in [[i for i in re.findall('[A-Z][^A-Z]*',x)] for x in torsions]])
        TypeList=[]
        for char in sitelist:
            if not TypeList.__contains__(char):
                TypeList.append(char)
            TypeCombos=[]
            for i,elem1 in enumerate(TypeList):
                for j,elem2 in enumerate(TypeList):
                    if i>=j:
                        TypeCombos.append([elem1,elem2])
        for i,j in TypeCombos:
            if ( (sitelist.count(i)+sitelist.count(j))/(1+ (i==j))   >2*(bonds.count(''.join([i,j]))+bonds.count(''.join([j,i])))):
                safelist.append(''.join([i,j]))
        return self.RemovePalindromes(set(safelist))
    
    def FindBondTriples(self,sitelist,bondlist):
        triples=[]
        for pair1 in bondlist:
            for pair2 in bondlist:
                central=set(pair1).intersection(set(pair2))
                if pair1 != pair2 and len(central)==1:
                    cid=list(central)[0]
                    if triples.count([pair2[pair2.index(cid)-1],cid, pair1[pair1.index(cid)-1]])==0:
                        triples.append([ pair1[pair1.index(cid)-1],cid,pair2[pair2.index(cid)-1] ])
        return triples
    
    def FindBondQuartets(self,sitelist,bondlist,triples):
        quartets=[]
        for pair in bondlist:
            for trip in triples:
                central1=set(pair).intersection(set(trip))
                #print pair,trip,central1
                if len(central1)==1:
                    central2=trip[1]
                    cid1=list(central1)[0]
                    cid2=central2
                    tripcentralid=trip.index(cid1)
                    if tripcentralid==0:
                        cid3=trip[2]
                        quartets.append(''.join([sitelist[i-1] for i in [ pair[pair.index(cid1)-1],cid1,cid2,cid3]]))
                    elif tripcentralid==2:
                        cid3=trip[0]
                        quartets.append(''.join([sitelist[i-1] for i in [ cid3,cid2,cid1,pair[pair.index(cid1)-1]]]))
        return self.RemovePalindromes(set(quartets))	
    
    def RemovePalindromes(self,r):
        r=list(r)
        for i,first in enumerate(r):
            for j,second in enumerate(r):
                #if j!=i and (first[::-1]==second or first==second):
                if j!=i and ([x for x in re.findall('[A-Z][^A-Z]*',first)]==[z for z in re.findall('[A-Z][^A-Z]*',second)][::-1] or first==second):
                    r.pop(j)
        return r

    def GetAAtoCGMap(self, aapointer, inchimap, peptide):
        SEQ=[i for i in peptide]
        maplist=[]
        offset=0
        for ii,aa in enumerate(SEQ):
            if aa == '*':
                offset=1
            elif aa != '*':
                val = [ [n+aapointer[ii]-1 for n in x] for x in inchimap[ii-offset]]
                for group in val:
                    maplist.append(group)
        return maplist	
    
    def GetCGInfo(self):
        bondsites=[]
        lastN=None
        siteid=[0]
        sitelist=[]
        inchimap=[]
        #remove end caps
        if not type(self.resseq)==type([]):
            peptide=peptide.replace('*','')
        else:
            if '*' in self.resseq: self.resseq.remove('*')
            if '*' in self.resseq: self.resseq.remove('*')
            SEQ=[i for i in self.resseq]
        for ii,aa in enumerate(SEQ):
            if not ii == 0:
                siteid.append(len(sitelist))
            inchimap.append(self.cg_dict[aa][3])
            inchi=self.cg_dict[aa][0].split('/')
            Nind= self.cg_dict[aa][1]
            Cind= self.cg_dict[aa][2]
            formula=inchi[1]
            bonds=inchi[2].lstrip('c')
            bondsites.extend(self.GetBondList(bonds,siteid))
            if not (ii in [x[0] for x in self.aapoint]):#ii==0:#^%$^%$
                bondsites.append([lastC+siteid[-2],Nind+siteid[-1]])
            lastC=Cind
            sitelist.extend(self.ParseInchiFormula(formula)[0])
        return bondsites,[i+1 for i in siteid],sitelist,inchimap
    
    
    def ParseInchiFormula(self,formula):
        charlist=[]
        charnumlist=[]
        for i,char in enumerate(formula):
            if char.isalpha() and not char.islower():
                StillLower=True
                while StillLower:
                    try:
                        if formula[i+1].islower():
                            char+=formula[i+1]
                            i=i+1
                        else:
                            StillLower=False
                    except: #exception occurs here if the final entry is an uppercase alpha character
                        pass
                        StillLower=False
                try:
                    if formula[i+1].isdigit():
                        NotAlpha=True
                        x=1
                        while NotAlpha:
                            try:
                                if formula[i+1+x].isalpha():
                                    NotAlpha=False
                                else:
                                    x+=1
                            except:
                                NotAlpha=False                   
                        charnumlist.append(int(formula[i+1:i+1+x]))
                        for j in range(int(formula[i+1:i+1+x])):
                            charlist.append(char)
                    else: #the next entry is a nondigit in which case it is capitol and a new entry 
                        charlist.append(char)
                        charnumlist.append(1)
                except:
                    #exception occurs if final formula entry is a alpha character
                    charlist.append(char)
                    charnumlist.append(1)
        return charlist,charnumlist
    
    def WriteMapFileString(self):
        fs="cgdict={"
        for k,v in self.cg_dict.items():
            fs+="""\n"%s" : %s,"""%(str(k),str(v))
        fs.rstrip(',')
        fs+="\n}\n\n"
        fs+="element_dict={"
        for k,v in self.element_dict.items():
            fs+="""\n"%s" : %s,"""%(str(k),str(v))
        fs.rstrip(',')
        fs+="\n}\n\n"
        return fs
    
    def showatom(self,molnum,num):
            self.Sys.Mol[molnum][num].Opacity=1.0
    
    def hideatom(self,molnum,num):
        self.Sys.Mol[molnum][num].Opacity=0.1
    
    def showmol(self,num):
        for atom in self.Sys.Mol[num]:
	        atom.Opacity=1.
    
    def showghostmol(self,num):
            for atom in self.Sys.Mol[num]:
                    atom.Opacity=.1
    
    def hidemol(self,num):
            for atom in self.Sys.Mol[num]:
                    atom.Opacity=.1
    
    def showall(self):
        for mol in self.Sys.Mol:
	        for atom in mol:
		        atom.Opacity=1.
    
    def hideall(self):
        for mol in self.Sys.Mol:
	        for atom in mol:
		        atom.Opacity=.0
    
    def deletesphere(self,x):
        self.Sys.Mol[0][x].Opacity=0.
    
    def InputGroups(self):
        Invalid=True
        display_error=0
        while Invalid:
            if display_error:
                print 'NOTE: CG group format is nested lists [[1,2],[3]]'
            s=str(raw_input('Input CG group atoms:\n'))
            try:
                l=eval(s)
                if not np.array(l).ndim==2:
                    display_error=1
                Invalid=False
            except:
                display_error=1
                pass
        return l
    
    def InputSiteName(self):
        Invalid=True
        while Invalid:
            try:
                s=str(raw_input('Sitename:\n'))
                s=s.upper()[0]+s.lower()[1:]
                if any(i.isdigit() for i in s):
                    for i in s:
                        print i,i.isdigit()
                    print 'NOTE: Sitenames cannot contain digits.'
                else:
                    Invalid=False
            except:
                pass
        return s   
    
    def GetProperties(self):
        if self.newsitename not in self.cgsite_dict:
            gcenter,gcolor,grg=self.GenProperties()
        else:
            gcenter,grg,gmass=self.GenCenter()
            grg,gcolor,mass=self.cgsite_dict[self.newsitename]
        return gcenter,gcolor,grg
    
    def GenProperties(self):
        gcenter,grg,gmass=self.GenCenter()
        gcolor=(random.random(),random.random(),random.random())
        self.cgsite_dict[self.newsitename]=[grg,gcolor,gmass]
        if self.newsitename == 'Wat':
            opacity=0.2
            gcolor=(0,0,1)
        else:
            opacity=1.0
        self.element_dict[self.newsitename]=[grg,gcolor,gmass,opacity]
        return gcenter,gcolor,grg
    
    def GenCenter(self):
        possum =np.zeros(self.Sys.Dim, dtype=float)
        possumsq=np.zeros(self.Sys.Dim, dtype=float)
        lastpos=np.zeros(self.Sys.Dim, dtype=float)
        masssum=0
        for atom in self.group:
            thispos=lastpos+self.Sys.MinImage(self.Sys.Mol[self.molid][self.reducedatoms[atom-1]].Pos-lastpos)
            possum+=thispos
            possumsq+=thispos**2.
            lastpos=thispos
            masssum+=self.element_dict[self.Sys.Mol[self.molid][self.reducedatoms[atom-1]].Name][2]
        gcenter=possum/len(self.group)
        grg=np.sqrt(np.sum(1./len(self.group)* possumsq - gcenter**2.))
        if grg == 0:
            grg=1.
            print 'default size',grg
        gmass=masssum
        return gcenter,grg,gmass
    
    def InputConnectivity(self):
        Invalid=True
        while Invalid:
            c=str(raw_input('Input connectivity of CG groups:\n'))
            try:
                thisbondlist=self.GetBondList(c)
                Invalid=False
            except:
                print 'Use InChI bond format.'
        for bonds in thisbondlist:
            self.cgrods.append(self.vis.v.cylinder(pos = self.cgspheres[bonds[0]-1].pos, axis = self.cgspheres[bonds[1]-1].pos-self.cgspheres[bonds[0]-1].pos, radius = .25/self.vscale, opacity = 0.5))
        Invalid=True
        while Invalid:
            inind=str(raw_input('Input id of head bead:\n'))
            try:
                int(inind)
                Invalid=False
            except:
                print 'Invalid.  Must be an interger'
        Invalid=True
        while Invalid:
            outind=str(raw_input('Input id of tail bead:\n'))
            try:
                int(outind)
                Invalid=False
            except:
                print 'Invalid.  Must be an interger'
        try:
            self.cg_dict[self.aa]=eval("['1S/"+self.formula+"/c"+c+"',"+inind+","+outind+","+str(self.groups)+"]")
            InvalidInchi=False
        except:
            print "['1S/"+self.formula+"/c"+c+"',"+inind+","+outind+","+str(self.groups)+"]"
            print 'Invalid Inchi please try again'
        for sphere in self.cgspheres:
            sphere.opacity=0.2
        for rod in self.cgrods:
            rod.opacity=0.2

    def GetBondList(self,bonds,siteid=None):
        bondsites=[]
        #find all hyphens and bond -1,+1...also...#find all '('  and bond -1,+1
        hyphens=[n for n in xrange(len(bonds)) if bonds.find('-', n) == n]
        lparens=[n for n in xrange(len(bonds)) if bonds.find('(', n) == n]
        hyphens.extend(lparens)
        for hyphen in hyphens:
            #find complete site numbers before the hyphen
            IsDigit=True
            x=1
            while IsDigit:
                if hyphen-x-1<0:
                    IsDigit=False
                else:
                    if bonds[hyphen-x-1].isdigit():
                        x+=1
                    else:
                        IsDigit=False
            before=bonds[hyphen-x:hyphen]
            #find complete site numbers after the hyphen
            IsDigit=True
            x=1
            while IsDigit:
                try:
                    if bonds[hyphen+x+1].isdigit():
                        x+=1
                    else:
                        IsDigit=False
                except:
                    IsDigit=False
            after=bonds[hyphen+1:hyphen+x+1]
            if siteid==None:
                bondsites.append([int(before),int(after)])
            else:
                bondsites.append([int(before)+siteid[-1],int(after)+siteid[-1]])
        #find all parenthesis pairs and for each ')' bond +1 and paired '('-1
        rparens=[n for n in xrange(len(bonds)) if bonds.find(')', n) == n]
        #depth=0
        indexlist=[]
        pairlist=[]
        for i in range(len(bonds)):
            if lparens.count(i)>=1: #if position is '('
                indexlist.append(i)
                #depth+=1
            if rparens.count(i)>=1:
                #depth-=1
                pairlist.append([indexlist[-1],i])
                indexlist.pop(-1)   
        for pair in pairlist: #for each pair or parenthesis find the bond targets
            #find complete site before the '('
            IsDigit=True
            x=1
            while IsDigit:
                if pair[0]-x-1<0:
                    IsDigit=False
                else:
                    if bonds[pair[0]-x-1].isdigit():
                        x+=1
                    else:
                        IsDigit=False
            before=bonds[pair[0]-x:pair[0]]
            #find complete site numbers after the paired ')'
            IsDigit=True
            x=1
            while IsDigit:
                try:
                    if bonds[pair[1]+x+1].isdigit():
                        x+=1
                    else:
                        IsDigit=False
                except:
                    IsDigit=False
            after=bonds[pair[1]+1:pair[1]+x+1]
            if siteid==None:
                bondsites.append([int(before),int(after)])
            else:
                bondsites.append([int(before)+siteid[-1],int(after)+siteid[-1]])
        #code to handle commas
        #find all '(' ',' pairs and for each ',' bond +1 and paired '('-1
        rparens=[n for n in xrange(len(bonds)) if bonds.find(')', n) == n]
        commas=[n for n in xrange(len(bonds)) if bonds.find(',', n) == n]
        #depth=0
        indexlist=[]
        commalist=[]
        for i in range(len(bonds)):
            if lparens.count(i)>=1:
                indexlist.append(i)
                #depth+=1
            if rparens.count(i)>=1:
                #depth-=1
                indexlist.pop(-1)   
            if commas.count(i)>=1:
                commalist.append([indexlist[-1],i])
        for comma in commalist:
            #find complete site before the '('
            IsDigit=True
            x=1
            while IsDigit:
                if comma[0]-x-1<0:
                    IsDigit=False
                else:
                    if bonds[comma[0]-x-1].isdigit():
                        x+=1
                    else:
                        IsDigit=False
            before=bonds[comma[0]-x:comma[0]]
            #find complete site numbers after the paired ')'
            IsDigit=True
            x=1
            while IsDigit:
                try:
                    if bonds[comma[1]+x+1].isdigit():
                        x+=1
                    else:
                        IsDigit=False
                except:
                    IsDigit=False
            after=bonds[comma[1]+1:comma[1]+x+1]
            if siteid==None:
                bondsites.append([int(before),int(after)])
            else:
                bondsites.append([int(before)+siteid[-1],int(after)+siteid[-1]])
        return bondsites

    
    def MakeMappingSys(self):
        framemap=[]
        for ii,i in enumerate(self.molecule):
            framemap.extend(i[-1])
        from sim import chem, traj, atomselect, potential, atommap, system
        #AtomTypes
        AtomTypes=[]
        for stype in self.sitetypes:
            if stype in self.element_dict:
                radius,color,mass,opacity=self.element_dict[stype]
            else:
                radius,color,mass,opacity=(np.random.random(),np.random.random(3),1,1)
            AtomTypes.append(chem.AtomType(stype,Mass=1., Charge=0., Color=(color[0],color[1],color[2]), Radius=radius, Opacity=opacity))
        
        #MolTypes
        MolTypes=[]
        nmol=len(self.networksitetypes)
        for ii,mol in enumerate(self.networksitetypes):
            MolTypes.append(chem.MolType( "molecule"+str(ii), [AtomTypes[self.site_dict[atomtype]] for atomtype in mol]  ))
        
        jj=ii
        for ii,mono in enumerate(self.monomers):
            MolTypes.append(chem.MolType( "molecule"+str(ii+jj+1), [AtomTypes[self.site_dict[mono]]]))
        
        #MolType Bonds
        for ii,bonds in enumerate(self.respooledbonds):
            minsite=min(min(bonds))
            for bond in bonds:
                MolTypes[ii].Bond(bond[0]-minsite,bond[1]-minsite)
        
        World=chem.World(MolTypes, Dim=3)
        SysName='Simulation'+str(time.time())[5:-2].rstrip('.')
        print 'creating system', SysName
        Sys = system.System(World, Name = SysName)
        
        if os.path.isfile(self.Prefix+'.mdtrj.crd.bz2') and os.path.isfile(self.Prefix+'.mdtrj.crd.bz2'):
            print 'loading all-atom trajectory'
            Trj=traj.Amber(self.Prefix+'.mdtrj.crd.bz2', self.Prefix+'.prmtop.parm7')
            frame=Trj[0]
            Sys.BoxL=Trj.FrameData['BoxL']
        else:
            print 'Error: no trajectory file found'
                
        for mol in MolTypes:
            Sys+=mol.New()
        
        PBond = potential.Bond(Sys, Filter = atomselect.BondPairs,
                       Dist0=0.8, FConst=3, Label = "Bond")
        Sys.ForceField.append(PBond)
        
        #set up the integrator 
        Sys.Int.Method = Sys.Int.Methods.VVIntegrate
        
        #define some things to measure.
        #setup a measurement of a histogram of potential energy
        Sys.Measures.PEnergy.SetupHist(-1200, -500, 300)
        
        #compile and load the system
        print 'loading system'
        Sys.Load()
        Sys.Pos=Trj[0][framemap]
        return Sys



        

class ReadGenerator(InterpretTopology):

    def __init__(self,Prefix):
        self.Prefix=Prefix
        namespace={}
        execfile(self.Prefix+'_generatorCG.txt',namespace)
        self.cgbondlist=namespace['cgbondlist']
        self.cgsitelist=namespace['cgsitelist']
        self.element_dict=namespace['element_dict']
        self.ambermap=namespace['ambermap']
        self.bondsets=namespace['bondsets']
        self.anglesets=namespace['anglesets']
        self.torsionsets=namespace['torsionsets']
        self.nonbondsets=namespace['nonbondsets']
        self.bondlist=self.cgbondlist
        self.sitelist=self.cgsitelist
        InterpretTopology.__init__(self)
        


    

    
    def MakeSysObject(self,Visualize=False,Rho=0.00001,TempSet=300.):
        from sim import chem, traj, atomselect, potential, atommap, system
        #AtomTypes
        AtomTypes=[]
        for stype in self.sitetypes:
            if stype in self.element_dict:
                radius,color,mass,opacity=self.element_dict[stype]
            else:
                radius,color,mass,opacity=(np.random.random(),np.random.random(3),1,1)
            AtomTypes.append(chem.AtomType(stype,Mass=1., Charge=0., Color=(color[0],color[1],color[2]), Radius=radius, Opacity=opacity))
        
        #MolTypes
        MolTypes=[]
        nmol=len(self.networksitetypes)
        for ii,mol in enumerate(self.networksitetypes):
            MolTypes.append(chem.MolType( "molecule"+str(ii), [AtomTypes[self.site_dict[atomtype]] for atomtype in mol]  ))
        
        jj=ii
        for ii,mono in enumerate(self.monomers):
            MolTypes.append(chem.MolType( "molecule"+str(ii+jj+1), [AtomTypes[self.site_dict[mono]]]))
        
        #MolType Bonds
        for ii,bonds in enumerate(self.respooledbonds):
            minsite=min(min(bonds))
            for bond in bonds:
                MolTypes[ii].Bond(bond[0]-minsite,bond[1]-minsite)
        
        #create atom filters
        FilterList=[]
        for ii,stype in enumerate(self.sitetypes):
            FilterList.append(atomselect.Filter(Types=AtomTypes[ii]))
        
        #create polyfilters
        BondFilters=[]
        for bondtype in self.bondsets:    
            #split the potentialtype string into useful components
            stypes=re.sub( r"([A-Z])", r" \1", bondtype).split()[1:]
            BondFilters.append(atomselect.PolyFilter(Filters=[ FilterList[self.site_dict[stype]] for stype in stypes],Bonded=True))
        
        AngleFilters=[]
        for angletype in self.anglesets:    
            #split the potentialtype string into useful components
            stypes=re.sub( r"([A-Z])", r" \1", angletype).split()[1:]
            AngleFilters.append(atomselect.PolyFilter(Filters=[ FilterList[self.site_dict[stype]] for stype in stypes],Bonded=True))
        
        TorsionFilters=[]
        for torsiontype in self.torsionsets:    
            #split the potentialtype string into useful components
            stypes=re.sub( r"([A-Z])", r" \1", torsiontype).split()[1:]
            TorsionFilters.append(atomselect.PolyFilter(Filters=[ FilterList[self.site_dict[stype]] for stype in stypes],Bonded=True))
        
        NonBondFilters=[]
        for nonbondtype in self.nonbondsets:    
            #split the potentialtype string into useful components
            stypes=re.sub( r"([A-Z])", r" \1", nonbondtype).split()[2:]
            NonBondFilters.append(atomselect.PolyFilter(Filters=[ FilterList[self.site_dict[stype]] for stype in stypes],Intra=True,MinBondOrd=4))
        
        World=chem.World(MolTypes, Dim=3)
        SysName='Simulation'+str(time.time())[5:-2].rstrip('.')
        print 'creating system', SysName
        Sys = system.System(World, Name = SysName)
        
        import os.path
        if os.path.isfile(self.Prefix+'.mdtrj.crd.bz2') and os.path.isfile(self.Prefix+'.mdtrj.crd.bz2'):
            print 'loading all-atom trajectory'
            Trj=traj.Amber(self.Prefix+'.mdtrj.crd.bz2', self.Prefix+'.prmtop.parm7')
            frame=Trj[0]
            Sys.BoxL=Trj.FrameData['BoxL']
        else:
            frame=None
            if Rho is not None:
                Sys.BoxL = (float(len(self.sitelist)) / Rho)**(1./Sys.Dim)
            else:
                Sys.BoxL = np.array([0.,0.,0.])
                
        for molid in self.molidlist:
            Sys+=MolTypes[molid].New()
        
        #create potentials
        BondPotentials=[]
        for ii,bondtype in enumerate(self.bondsets):    
            stypes=re.sub( r"([A-Z])", r" \1", bondtype).split()[1:]
            BondPotentials.append(potential.Bond(Sys, Label = "%s"%bondtype, Filter = BondFilters[ii], Dist0= 4.0, FConst = 1.0))
        
        AnglePotentials=[]
        for ii,angletype in enumerate(self.anglesets):    
            stypes=re.sub( r"([A-Z])", r" \1", angletype).split()[1:]
            AnglePotentials.append(potential.AngleSpline(Sys, Label = "%s"%angletype, Filter = AngleFilters[ii], NKnot=20))
        
        TorsionPotentials=[]
        for ii,torsiontype in enumerate(self.torsionsets):    
            stypes=re.sub( r"([A-Z])", r" \1", torsiontype).split()[1:]
            TorsionPotentials.append(potential.TorsionSpline(Sys, Label = "%s"%torsiontype, Filter = TorsionFilters[ii], NKnot=20))
        
        NonBondPotentials=[]
        for ii,nonbondtype in enumerate(self.nonbondsets):    
            stypes=re.sub( r"([A-Z])", r" \1", nonbondtype).split()[1:]
            NonBondPotentials.append(potential.PairSpline(Sys, Label = "%s"%nonbondtype, Cut = 15.0, Filter = NonBondFilters[ii], NKnot=40))
        
        Sys.ForceField.extend([j for i in [BondPotentials,AnglePotentials,TorsionPotentials,NonBondPotentials] for j in i])
        #set up the integrator 
        Sys.Int.Method = Sys.Int.Methods.VVIntegrate
        
        #define some things to measure.
        #setup a measurement of a histogram of potential energy
        Sys.Measures.PEnergy.SetupHist(-1200, -500, 300)
        
        #compile and load the system
        print 'loading system'
        Sys.Load()
        
        #create aa->cg mapping rule
        Sys.Map = atommap.PosMap()
        i=0 
        for ii,mol in enumerate(self.molidlist):
            if mol>=nmol:
                natom=1
            else:
                natom=len(self.networksitetypes[mol])
            for jj in range(natom):
                Sys.Map += [atommap.AtomMap( Atoms1 = [k-1 for k in self.ambermap[i]], Atom2 = Sys.Mol[ii][jj] )] 
                i+=1
        
        if frame is None:
            if Rho is not None:
                print 'setting positions to cubic lattice'
                system.init.CubicLattice(Sys, Random = 0.1)
            else:
                print 'setting positions to origin: Note that positions need to be set manually'
        else:
            print 'setting positions to first trajectory frame'
            MappedTrj=traj.Mapped(Trj,Sys.Map)
            Sys.Pos=MappedTrj[0]
        
        system.velocities.Canonical(Sys, Temp = TempSet)

        if Visualize:
            import sim.system.visualize as vis
            Sys.Vis = vis.Visualizer3D(Sys,ShowBonds=True,Label=0)
            Sys.Vis.Update()
            Sys.Vis.AddAction(Sys.Int, TimeFreq = 0.2)
            Sys.Measures.VerboseOutput(StepFreq = 1000)


        print 'initializeing verlet integrator'    
        Sys.Int.Method = Sys.Int.Methods.VVIntegrate
        return Sys
    
    def MakeSysScript(self, Visualize=False, Rho=0.00001, TempSet=300.):
        ms="""
from sim import *
element_dict={"""
        for k,v in self.element_dict.items():
            ms+="""\n"%s" : %s,"""%(str(k),str(v))
        ms.rstrip(',')
        ms+="\n}"
        #AtomTypes
        ms+="\n\nprint 'building atom types'\nAtomTypes=[ "
        AtomTypes=[]
        for stype in self.sitetypes:
            if stype in self.element_dict:
                radius,color,mass,opacity=self.element_dict[stype]
            else:
                radius,color,mass,opacity=(np.random.random(),np.random.random(3),1,1)
            #AtomTypes.append(chem.AtomType(stype,Mass=1., Charge=0., Color=color, Radius=radius))
            ms+="""\nchem.AtomType("%s",Mass=1., Charge=0., Color=(%s,%s,%s), Radius=%s, Opacity=%s),"""%(stype,color[0],color[1],color[2],radius,opacity)
        
        ms=ms.rstrip(',')
        ms+="]\n\n"
        
        #MolTypes
        ms+="\nprint 'building moltypes'\nMolTypes=[ "
        MolTypes=[]
        nmol=len(self.networksitetypes)
        for ii,mol in enumerate(self.networksitetypes):
            ms+="""\nchem.MolType( "molecule%d", ["""%(ii)
            for atomtype in mol:
                ms+=""" AtomTypes[%d],"""%self.site_dict[atomtype]
            
            ms=ms.rstrip(",")
            ms+="]),"
        
        jj=ii
        for ii,mono in enumerate(self.monomers):
            ms+="""\nchem.MolType( "molecule%d", [AtomTypes[%d]]),"""%(ii+jj+1,self.site_dict[mono] )
        
        ms=ms.rstrip(",")
        ms+="]\n"
        
        
        ms+="\nprint 'bonding moltypes'"
        #MolType Bonds
        for ii,bonds in enumerate(self.respooledbonds):
            minsite=min(min(bonds))
            for bond in bonds:
                ms+="\nMolTypes[%d].Bond(%d,%d)"%(ii,bond[0]-minsite,bond[1]-minsite)
        
        #create atom filters
        ms+="""\n\nprint "creating atom filters"\nFilterList=["""
        for ii,stype in enumerate(self.sitetypes):
            ms+="\natomselect.Filter(Types=AtomTypes[%d]),"%ii
        
        ms=ms.rstrip(",")
        ms+="]\n\n"
        
        #create polyfilters
        #bondfilters
        ms+="BondFilters=["
        for bondtype in self.bondsets:    
            #split the potentialtype string into useful components
            stypes=re.sub( r"([A-Z])", r" \1", bondtype).split()[1:]
            ms+="\natomselect.PolyFilter(Filters=["
            for stype in stypes:
                ms+="FilterList[%d],"%self.site_dict[stype]
            
            ms=ms.rstrip(",")
            ms+="],Bonded=True),"
        
        ms=ms.rstrip(",")
        ms+="]\n\n"
        
        ms+="AngleFilters=["
        for angletype in self.anglesets:    
            #split the potentialtype string into useful components
            stypes=re.sub( r"([A-Z])", r" \1", angletype).split()[1:]
            ms+="\natomselect.PolyFilter(Filters=["
            for stype in stypes:
                ms+="FilterList[%d],"%self.site_dict[stype]
            
            ms=ms.rstrip(",")
            ms+="],Bonded=True),"
        
        ms=ms.rstrip(",")
        ms+="]\n\n"
        
        ms+="TorsionFilters=["
        for torsiontype in self.torsionsets:    
            #split the potentialtype string into useful components
            stypes=re.sub( r"([A-Z])", r" \1", torsiontype).split()[1:]
            ms+="\natomselect.PolyFilter(Filters=["
            for stype in stypes:
                ms+="FilterList[%d],"%self.site_dict[stype]
            
            ms=ms.rstrip(",")
            ms+="],Bonded=True),"
        
        ms=ms.rstrip(",")
        ms+="]\n\n"
        
        ms+="NonBondFilters=["
        for nonbondtype in self.nonbondsets:    
            #split the potentialtype string into useful components
            stypes=re.sub( r"([A-Z])", r" \1", nonbondtype).split()[2:]
            ms+="\natomselect.PolyFilter(Filters=["
            for stype in stypes:
                ms+="FilterList[%d],"%self.site_dict[stype]
            
            ms=ms.rstrip(",")
            ms+="], Intra = True, MinBondOrd = 4),"
        
        ms=ms.rstrip(",")
        ms+="]\n\n"
        
        #Start building the system object
        ms+="""
print 'creating system'
World=chem.World(MolTypes, Dim=3)
SysName='Simulation'
Rho=%11.5f
TempSet=%11.5f
Sys = system.System(World, Name = SysName)
"""%(Rho,TempSet)
        
        
        #load aa traj if possible
        ms+="""

print 'loading all-atom trajectory'
Trj=traj.Amber('%s.mdtrj.crd.bz2', '%s.prmtop.parm7')
frame=Trj[0]
Sys.BoxL=Trj.FrameData['BoxL']
        """%(self.Prefix,self.Prefix)
        
        #add in molecules 
        ms+="""

print 'adding molecules to system'
molidlist=%s
for molid in molidlist:
    Sys+=MolTypes[molid].New()

"""%str(self.molidlist)
        
        #ms+="""
        #set the system box length sizes
        #Sys.BoxL = (float(len(molidlist)) / Rho)**(1./Sys.Dim)
        #"""
        #create potentials
        ms+="BondPotentials=["
        for ii,bondtype in enumerate(self.bondsets):    
            stypes=re.sub( r"([A-Z])", r" \1", bondtype).split()[1:]
            ms+="""\npotential.Bond(Sys, Label = "%s", Filter = BondFilters[%d], Dist0= 4.0, FConst = 1.0),"""%(bondtype,ii)
        
        ms=ms.rstrip(",")
        ms+="]\n\n"
        
        ms+="AnglePotentials=["
        for ii,angletype in enumerate(self.anglesets):    
            stypes=re.sub( r"([A-Z])", r" \1", angletype).split()[1:]
            ms+="""\npotential.AngleSpline(Sys, Label = "%s", Filter = AngleFilters[%d], NKnot=20),"""%(angletype,ii)
        
        ms=ms.rstrip(",")
        ms+="]\n\n"
        
        ms+="TorsionPotentials=["
        for ii,torsiontype in enumerate(self.torsionsets):    
            stypes=re.sub( r"([A-Z])", r" \1", torsiontype).split()[1:]
            ms+="""\npotential.TorsionSpline(Sys, Label = "%s", Filter = TorsionFilters[%d], NKnot=20),"""%(torsiontype,ii)
        
        ms=ms.rstrip(",")
        ms+="]\n\n"
        
        ms+="NonBondPotentials=["
        for ii,nonbondtype in enumerate(self.nonbondsets):    
            stypes=re.sub( r"([A-Z])", r" \1", nonbondtype).split()[1:]
            ms+="""\npotential.PairSpline(Sys, Label = "%s", Cut = 15.0, Filter = NonBondFilters[%d], NKnot=40),"""%(nonbondtype,ii)
        
        ms=ms.rstrip(",")
        ms+="]\n\n"
        
        ms+="""
Sys.ForceField.extend([j for i in [BondPotentials,AnglePotentials,TorsionPotentials,NonBondPotentials] for j in i])
        """
        
        
        
        ms+="""

#set up the integrator 
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate

#define some things to measure.
#setup a measurement of a histogram of potential energy
Sys.Measures.PEnergy.SetupHist(-1200, -500, 300)

#compile and load the system
print 'loading system'
Sys.Load()
        """
        ms+="""
Sys.Map = atommap.PosMap()"""
        
        
        i=0 
        for ii,mol in enumerate(self.molidlist):
            if mol>=nmol:
                natom=1
            else:
                natom=len(self.networksitetypes[mol])
            for jj in range(natom):
                ms+="\nSys.Map += [atommap.AtomMap( Atoms1 = %s, Atom2 = Sys.Mol[%d][%d] )]"%(str([k-1 for k in self.ambermap[i]]), ii, jj) 
                i+=1
        
        
        ms+="""
MappedTrj=traj.Mapped(Trj,Sys.Map)
Sys.Pos=MappedTrj[0]

import sim.system.visualize as vis
Sys.Vis = vis.Visualizer3D(Sys,ShowBonds=True,Label=0)
Sys.Vis.Update()
Sys.Vis.AddAction(Int, TimeFreq = 0.2)
Sys.Measures.VerboseOutput(StepFreq = 1000)

TempSet = 300.0 #Kelvin
system.velocities.Canonical(Sys, Temp = TempSet)
#Int.Method = Int.Methods.VVQuench
#print "Minimizing"
#Int.Run(2000)

#Int.Method = Int.Methods.VVIntegrate
#Sys.TempSet = TempSet
#Int.Run(2000)
        """
        
        fo =open(self.Prefix+'_md.py','w')
        fo.write(ms)
        fo.close()
        
        return
        































