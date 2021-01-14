#!/usr/bin/env python
import glob
import numpy as np
import scipy.interpolate as fit

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


DEBUG = False

#this is set off by default
SmoothHistData = False

#outlier handling
OutlierZScore = 10.0     

#maximum labels or won't put a legend
LegendMaxLabels = 5

#plots per row for iterations
IterationsPlotsPerRow = 5

#max nonbonded ene
MaxNonbondEne = "5kT"
MaxParamVal = 5.

#show all type combinations?
ReportAllTypes = True


class PlotForceField(object):
    
    def __init__(self, Sys, filename = None, modhistfile = None, modhistind = 0):
        self.Sys=Sys
        self.fig=None
        self.ax=None
        self.xdim=None
        self.ydim=None
        self.spacevar = 1
        self.xadjust = 1
        self.pot_x=None
        self.pot_y=None
        self.fluct_y = None
        self.tar_x=None
        self.tar_y=None
        self.mod_x=None
        self.mod_y=None
        self.MaxNonbondEne = Sys.Units.ConvertEne(MaxNonbondEne, Sys.TempSet)
        self.GeneratePlot(filename = filename, modhistfile = modhistfile, modhistind = modhistind)
        return
        
    def getFile(self, ftype):
        #fetch a file and return it as a file object
        fl=glob.glob(ftype)
        if len(fl)>1 or len(fl)==0:
            raise IndexError("Too many files, or no files")
        return open(fl[0],'r')
        
    def parseModHist(self, filename='*modhist.txt'):
        #opens the modhist txt file and parses out the histogram information
        #returns a dictionary with key=potential_name and val= Nx3 array of x,target,model
        fobj = self.getFile(filename)
        fread = fobj.read().strip()
        histograms = fread.split('POTENTIAL ')
        mod_dict = {}
        for block in histograms[1:]:
            data = []
            lines = block.strip().split("\n")
            for line in lines[3:]:
                data.append(np.fromstring(line, dtype=float, sep=' '))
            try:
                data = np.array(data, dtype=float)
            except ValueError:
                print "Error parsing mod hist file:"
                print data
                raise ValueError
            #split data into columns
            labels = lines[2].lower().strip().split()
            data = list(data.T)
            newdata = []
            #add in missing argument information
            lastarg = data[0]
            lastargused = False
            for i in range(len(data)):
                if labels[i] == 'arg': 
                    lastarg = data[i]
                    lastargused = False
                elif lastargused:
                    #need to add this arg column again
                    newdata.append(lastarg)
                else:
                    lastargused = True
                newdata.append(data[i]) 
            newdata = np.array(newdata).T
            name = lines[0].strip()
            mod_dict[name] = newdata
        return mod_dict
        
    def makeAxisLabel(self, pot, name):
        import sim.potential.base.potentialtypes as ptypes
        #this will return a nicely formated axis label with units that match the potential type
        if pot.Type == ptypes.PairPotential:
            unit = self.Sys.Units.LLabel
            if len(unit) > 0: unit = " (" + unit.strip() + ")"
            if pot.Filter.Bonded:
                return name + " bond dist" + unit
            else:
                return name + " pair dist" + unit
        elif pot.Type == ptypes.AnglePotential:
            return name + r" bond angle $(\theta)$"
        elif pot.Type == ptypes.TorsionPotential:
            return name + r" torsion angle $(\psi)$"
        else:
            return name + " arg" 
        
    def fixMarkers(self,ax):
        #resizes the tick marks and formats their labels
        ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(18)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(18)
        for l in ax.get_xticklines():
            l.set_markersize(6)
        for l in ax.get_yticklines():
            l.set_markersize(6)
        return
        
    def subplotDim(self):
        #this is an attempt to layout the subplots in some logical manner, but it needs improvement
        if self.dim<=4:
            self.xdim=2
            self.ydim=2
        elif self.dim==10:
            self.xdim=3
            self.ydim=4
        elif self.dim>=11 and self.dim <= 24:
            self.xdim=4
            self.ydim=6
        elif self.dim>24 and self.dim <=30:
            self.xdim=5
            self.ydim=6
        else:
            self.xdim=round(self.dim**.5)+1
            self.ydim=round(self.dim**.5)+1
        return
              
    def setRange(self, pot_element, typeind = 0):
        import sim.potential.base.potentialtypes as ptypes
        #check for histogram data
        HasHist = not getattr(self, "tar_x", None) is None
        #set x range
        if HasHist:
            self.xmin = self.tar_x.min()
            self.xmax = self.tar_x.max()
        else:
            self.xmin = pot_element.Arg.ReportMin
            self.xmax = pot_element.Arg.ReportMax
        #make potential values
        pot_element.SetTypeInd(typeind)
        self.pot_x = np.linspace(self.xmin, self.xmax, 100)
        self.pot_y = np.array([pot_element.Val(i) for i in self.pot_x])
        if pot_element.UsesFluct:
            self.fluct_y = np.array([pot_element.FLVal(i) for i in self.pot_x])
            self.fluct_y = np.clip(self.fluct_y, 0., None)
            self.fluct_y = np.sqrt(self.fluct_y)
        #make y range
        self.ymin = self.pot_y.min()
        self.ymax = self.pot_y.max()
        #check for nonbond potentials
        if pot_element.Type == ptypes.PairPotential:
            if not pot_element.Filter.Bonded:
                self.ymax = min(self.ymax, self.MaxNonbondEne)   
        ybuffer = abs(self.ymax-self.ymin) * 0.01
        self.ymin -= ybuffer
        self.ymax += ybuffer                
        self.ax.set_xlim(self.xmin, self.xmax)
        self.ax.set_ylim(self.ymin, self.ymax)
        self.ax.xaxis.set_major_locator(MaxNLocator(4))
        if HasHist:
            #the masks tell which indices in the data arrays are within xmin and xmax
            tar_mask = np.logical_and(self.tar_x <= self.xmax, self.tar_x >= self.xmin)
            mod_mask = np.logical_and(self.mod_x <= self.xmax, self.mod_x >= self.xmin)
            if DEBUG:
                print pot_element.Name
                print self.xmin, self.xmax
                n = len(self.tar_x)
                for i in range(n):
                    print "%5d %11.4e %11.4e %1d %11.4e %11.4e %1d" % (i, self.tar_x[i], self.tar_y[i], tar_mask[i], 
                                                                       self.mod_x[i], self.mod_y[i], mod_mask[i])
            if not(np.any(mod_mask)):
                self.ymin2 = self.tar_y[tar_mask,1].min()
                self.ymax2 = self.tar_y[tar_mask,1].max()            
            else:        
                self.ymin2 = min(self.tar_y[tar_mask].min(), self.mod_y[mod_mask].min())
                self.ymax2 = max(self.tar_y[tar_mask].max(), self.mod_y[mod_mask].max())
            ybuffer2 = abs(self.ymax2 - self.ymin2) * 0.01
            self.ymin2 = 0.
            self.ymax2 += ybuffer2
            self.ax2.set_ylim(self.ymin2, self.ymax2)           
        return self.pot_x, self.pot_y
               
    def GeneratePlot(self, filename, modhistfile=None, modhistind=0):
        if modhistfile:
            histograms = self.parseModHist(modhistfile)
        #make a list of forcefield elements
        Potentials = self.Sys.ForceField.GetUniquePotentials(AllTypes = ReportAllTypes)
        Potentials = sorted(Potentials, key = lambda x :x[1])
        self.dim = len(Potentials)
        #size plots
        self.subplotDim()
        #initialize figure and start drawing plots
        self.fig = plt.figure(figsize=(((self.xdim+self.xadjust)*5, self.ydim*5)))
        for (j, (P, PName, typeind)) in enumerate(Potentials):
            #make a plot
            self.ax = plt.subplot(self.ydim, self.xdim, j+1)
            if modhistfile:
                self.ax2 = self.ax.twinx()
                if modhistind == -1:
                    #average the histograms -- assums all have same arg bins
                    self.tar_x = histograms[PName][:,0::4].mean(axis=1)
                    self.tar_y = histograms[PName][:,1::4].mean(axis=1)
                    self.mod_x = histograms[PName][:,2::4].mean(axis=1)
                    self.mod_y = histograms[PName][:,3::4].mean(axis=1)
                else:
                    ind = modhistind*4
                    self.tar_x = histograms[PName][:,ind]
                    self.tar_y = histograms[PName][:,ind+1]
                    self.mod_x = histograms[PName][:,ind+2]
                    self.mod_y = histograms[PName][:,ind+3]
            #set the range of the axes
            self.pot_x, self.pot_y = self.setRange(P, typeind = typeind)
            #plot the potentials and label the axes
            p1 = self.ax.plot(self.pot_x, self.pot_y, label='potential', linewidth=3.0, color='red')
            #TODO only label subplots that fall along the outer edges to conserve space 
            self.ax.set_xlabel(self.makeAxisLabel(P, PName), size='xx-large')
            ylabel = 'potential energy (%s)' % self.Sys.Units.ELabel
            self.ax.set_ylabel(ylabel, size='xx-large', color='red')
            self.fixMarkers(self.ax)
            if P.UsesFluct:
                p2 = self.ax.plot(self.pot_x, self.fluct_y, label='fluctuations', linewidth=2.0, color='orange')
            #plot histograms
            if modhistfile:
                #plot the forcefield and histograms over this range
                p3 = self.ax2.plot(self.tar_x, self.tar_y, label='all atom', linewidth=3., color='black')
                p4 = self.ax2.plot(self.mod_x, self.mod_y, label='coarse-grained', linewidth=5., linestyle='--', color='blue')
                #label the axis
                self.ax2.set_ylabel('distribution', size="xx-large", color='blue')
                self.fixMarkers(self.ax2)
                #put a legend on the first subplot
                if j == 0:
                    self.leg = plt.legend() 
                    self.leg.get_frame().set_linewidth(0)
                    self.leg.get_frame().set_alpha(0.75)
        #final adjustments
        plt.subplots_adjust(wspace=0.7*self.spacevar, hspace=0.4*self.spacevar)
        try:
            plt.savefig(filename, bbox_inches='tight')
        except IOError:
            print "WARNING: Could not save histograms and potentials figure."
        plt.close()
        return
        


def CalculateYRanges(Data, yMin = 1.e300, yMax = -1.e300):
    Data = np.array(Data)
    while len(Data) > 5:
        i = np.argmax(np.abs(Data - Data.mean()))
        NewData = np.delete(Data, i)
        ZScore = np.abs(Data[i] - NewData.mean()) / NewData.std()
        if ZScore > OutlierZScore:
            Data = NewData
        else:
            yMin = min(yMin, Data.min())
            yMax = max(yMax, Data.max())
            yMax = min(yMax, MaxParamVal)
            return yMin, yMax
    yMax = min(yMax, MaxParamVal)
    return yMin, yMax
        

def PlotParameterHistory(logfile, filename, Skipdx = True):
    #make a dictionary of potentials; each will contain another
    #dictionary of parameters as well as "ITERS" which is a list
    d = {}
    #list of minimization starts
    StartIters = []
    IsStartIter = False
    #parse file
    LastIter = None
    for line in file(logfile, "r"):
        vals = line.strip().split()
        if len(vals) <=1:
            continue
        if line.lower().startswith("iter"):
            #this is a header line
            IsStartIter = True  
            #find name of potential and parameter label
            Params = []
            for Header in vals[2:]:
                if ":" in Header:
                    Potential, Label = Header.split(":")
                else:
                    Potential, Label = "Optimization", Header
                if any([x == (Potential, Label) for x in Params]):
                    raise ValueError("Found duplicate potential name %s when trying to plot forcefield." % Potential)
                Params.append((Potential, Label))
        else:
            #add more items to the parameter dictionaries
            Iter = int(vals[0])
            if IsStartIter:
                StartIters.append(Iter)
                IsStartIter = False
            elif vals[1] == "F":
                #flag for next minimization
                IsStartIter = True
            #check to see if this is a repeat
            if Iter == LastIter:
                continue
            LastIter = Iter
            #trim off values
            vals = vals[2:]
            for ((Potential, Label), val) in zip(Params, vals):
                #check if this is a constrained parameter; if so, skip
                if "*" in Label:
                    continue
                if Label == "dx" and Skipdx:
                    continue
                #convert parameter value to float
                val = float(val)
                #check if potential is yet defined
                if Potential in d:
                    Iters = d[Potential]["ITERS"]
                    if not Iters[-1] == Iter:
                        Iters.append(Iter)
                else:
                    d[Potential] = {"ITERS":[Iter]}
                #check if parameter is yet defined
                if Label in d[Potential]:
                    d[Potential][Label].append(val)
                else:
                    d[Potential][Label] = [val]
    #check that something was added to the dictionary
    if not len(d):
        raise ValueError("Nothing found to plot in log file %s" % logfile)
    #convert to numpy arrays
    for (Potential, dPotential) in d.iteritems():
        for Label in dPotential:
            dPotential[Label] = np.array(dPotential[Label], dtype=float)
    #set up the plotting
    N = len(d) 
    xdim = min(N, IterationsPlotsPerRow)
    ydim = int(np.ceil(float(N) / IterationsPlotsPerRow))
    fig = plt.figure(figsize=(xdim*5, ydim*5)) 
    #loop through potentials
    Potentials = ['Optimization'] + sorted([p for p in d.keys() if not p=='Optimization'])
    for (i, Potential) in enumerate(Potentials):
        sp = plt.subplot(ydim, xdim, i+1)
        sp.set_title(Potential)
        thisd = d[Potential]
        Iters = thisd.pop("ITERS")
        Labels = sorted(thisd.keys())
        yMin, yMax = 1.e300, -1.e300
        for Label in Labels:
            Vals = thisd[Label]
            yMin, yMax = CalculateYRanges(Vals, yMin, yMax)
            sp.plot(Iters, Vals, label=Label)
        if i == 0:
            sp.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        sp.set_xlim(0, LastIter)
        sp.set_ylim(yMin, yMax)
        if len(Labels) <= LegendMaxLabels:
            leg = plt.legend()
            leg.get_frame().set_linewidth(0)
            leg.get_frame().set_alpha(0.75)
        for StartIter in StartIters[1:]:
            plt.axvline(x=StartIter, color='k', linestyle=":")
    #do some finishing
    plt.subplots_adjust(wspace = 0.2, hspace = 0.2)
    #save
    try:
        plt.savefig(filename, bbox_inches='tight')
    except IOError:
        print "WARNING: Could not save parameter history figure."
    plt.close()

            
  