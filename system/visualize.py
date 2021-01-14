#!/usr/bin/env python

#Simple routines for visualizing a simulation of atoms.
#Requires the VPython module (version >= 5.0).

try:
    import visual as v
    HASVISUAL = True
except ImportError:
    print "Visual Python not available for molecular visualization."
    HASVISUAL = False

import numpy as np
import atexit



if HASVISUAL:
    
    class Visualizer3D(object):
    
        def __init__(self, Sys, BoxColor = v.color.blue, BackColor = v.color.white,
                     Scale = None, CloseOnExit = True, ShowBonds = False, 
                     Stereo = False, Label = False,
                     BondColor = (0.39, 0.39, 0.39), BondRadius = 0.25,
                     BoxRadius = 0.0001):
            """Initializes a 3D visualization window.
        Input:
            Sys: System object
            BoxColor: (R,G,B) tuple color for box; use None for no box
            BackColor: (R,G,B) tuple color for background
            Scale: length scale to which to normalize
            CloseOnExit: True to close visual window at program exit
            ShowBonds: Turn bond visualization on or off
            Stereo: Use stereo visualization
            Label: Label atoms with indices
            BondColor: Color of bonds
            BondRadius: relative size of bond radii
            BoxRadius: relative size of box boundary lines
        """
    
            #check dim
            if not Sys.Dim == 3:
                raise TypeError("System must have 3 dimensions")
    
            self.Sys = Sys        
            
            #save box length and scale
            self.__LastBoxL = Sys.BoxL.copy()
            self.__UseL = np.all(Sys.BoxL > 0)
            if Scale is None:
                if self.__UseL:
                    self.Scale = np.max(Sys.BoxL)
                else:
                    self.Scale = Sys.Units.LScale
            else:
                self.Scale = Scale
            self.iScale = 1. / self.Scale
            
            #minimum image and scale positions
            ScPos = self.Sys.Pos.copy()
            if self.__UseL:
                ScPos = ScPos - self.Sys.BoxL * np.around(ScPos / self.Sys.BoxL)
            ScPos = ScPos * self.iScale
    
            #lists of objects
            self.RodList = []
            self.BoxList = []
            self.SphereList = []
            self.Label = Label
            self.LabelList = []
            if (Stereo):
                self.Display = v.display(stereo = "redcyan")
            else:
                self.Display = v.display()
            self.Display.exit = False
            self.ShowBonds = ShowBonds
            
            #save colors and sizes
            self.BoxColor = BoxColor
            self.BoxRadius = BoxRadius
            self.BackColor = BackColor
            self.BondColor = BondColor
            self.BondRadius = BondRadius
           
            #decide whether or not to use the box
            self.__UseBox = self.__UseL and not BoxColor is None
            
            #draw everything
            self.Update()
    
            #change the background color
            self.Display.background = BackColor
    
            #set autocentering and scaling
            if not self.__UseBox:
                self.Display.autocenter = 1
                self.Display.autoscale = 1
                
            #check for close on exit
            if CloseOnExit:
                atexit.register(self.Close)
    
    
        def Update(self, Sys = None):
            """Updates a 3D visualization window.
Sys is here because Integrator actions will send a sys object."""
            self.Display.select()
            
            #check that sys is right
            if Sys is None:
                Sys = self.Sys
            elif not Sys is self.Sys:
                raise TypeError("Incorrect system object sent to Update.")
            
            #see if we need to update box
            if self.__UseBox:
                #check if we need to make the box
                #make a list of all corner points
                if len(self.BoxList) == 0:
                    Corners = []
                    for x in [-0.5*Sys.BoxL[0], 0.5*Sys.BoxL[0]]:
                        for y in [-0.5*Sys.BoxL[1], 0.5*Sys.BoxL[1]]:
                            for z in [-0.5*Sys.BoxL[2], 0.5*Sys.BoxL[2]]:
                                Corners.append((x,y,z))
                    Corners = np.array(Corners, float)
                    #sort through pairs of corners; for unit length, draw a cylinder
                    Radius = np.min(Sys.BoxL) * self.BoxRadius
                    for (i, c1) in enumerate(Corners):
                        for c2 in Corners[i+1:]:
                            if np.sum(c1 == c2) < 2: continue
                            self.BoxList.append(v.cylinder(pos = c1 * self.iScale,
                                                           axis = (c2-c1) * self.iScale,
                                                           radius = Radius,
                                                           color = self.BoxColor, opacity=1.0))
                else:
                    #update old box
                    if not np.all(self.__LastBoxL == Sys.BoxL):
                        Scale = Sys.BoxL / self.__LastBoxL
                        for BoxCyl in self.BoxList:
                            BoxCyl.axis = (Scale * BoxCyl.axis)
                            BoxCyl.pos = (Scale * BoxCyl.pos)
                #save the last box size
                self.__LastBoxL = Sys.BoxL.copy()
                    
            #minimum image and scale positions
            ScPos = Sys.Pos.copy()
            if self.__UseL:
                ScPos = ScPos - Sys.BoxL * np.around(ScPos / Sys.BoxL)
            ScPos = ScPos * self.iScale
    
            #update bonds
            if self.ShowBonds:
                RodInd = 0
                NRod = len(self.RodList)
                for (MolInd, Mol) in enumerate(Sys.Mol):
                    StartInd = Sys.MolRange[MolInd]
                    for Bond in Mol.MType.Bonds:
                        Ind1 = StartInd + Bond.SType1.AInd
                        Ind2 = StartInd + Bond.SType2.AInd
                        Pos1 = ScPos[Ind1]
                        Pos2 = ScPos[Ind2]
                        Atom1 = Sys.Atom[Ind1]
                        #minimium images rod axis
                        if self.__UseL:
                            Axis = (Pos2 - Pos1) - Sys.BoxL*self.iScale * np.around((Pos2 - Pos1) / (Sys.BoxL*self.iScale))
                        else:
                            Axis = (Pos2 - Pos1)
                        #check if this rod exists already
                        if RodInd < NRod:
                            #update existing rod
                            self.RodList[RodInd].pos = Pos1
                            self.RodList[RodInd].axis = Axis
                            self.RodList[RodInd].opacity = Atom1.Opacity
                        else:
                            #make a new rod
                            self.RodList.append(v.cylinder(pos = Pos1, axis = Axis, 
                                                           radius = self.BondRadius * self.iScale,
                                                           color = self.BondColor, opacity=Atom1.Opacity))                
                        #update counter
                        RodInd += 1
                
                #see if we need to delete any rods in case their number changed
                for i in range(RodInd, NRod):
                    self.RodList[-1].visible = False
                    del self.RodList[-1]
    
            #update positions of existing spheres
            NSphere = len(self.SphereList)
            for (i, p) in enumerate(ScPos):
                if i < NSphere:
                    #update existing sphere
                    self.SphereList[i].pos = p
                    self.SphereList[i].opacity = Sys.Atom[i].Opacity
                    if self.Label:
                        self.LabelList[i].pos = p
                else:
                    #make a new sphere
                    self.SphereList.append(v.sphere(pos = p,
                                           radius = Sys.Atom[i].Radius * self.iScale,
                                           color = Sys.Atom[i].Color,
                                           opacity = Sys.Atom[i].Opacity))
                    if self.Label:
                        self.LabelList.append(v.label(pos = p, text=str(Sys.Atom[i].Ind+1), box=0))
    
            #see if we need to delete any spheres in case their number changed
            for i in range(len(ScPos), NSphere):
                self.SphereList[-1].visible = False
                del self.SphereList[-1]
                if self.Label:
                    del self.LabelList[-1]
   
  
        def Clear(self):
            """Clears objects from the visualization window."""
            while len(self.SphereList):
                self.SphereList.pop(-1).visible = False
            while len(self.BoxList):
                self.BoxList.pop(-1).visible = False
            while len(self.RodList):
                self.RodList.pop(-1).visible = False
    
        def Close(self):
            """Closes the 3D visualization window."""
            self.Clear()
            if not self.Display is None:
                self.Display.exit = False
                self.Display.visible = False
                self.Display = None
    
        def AddAction(self, IntegrateObj, TimeFreq = 0.1):
            """Adds an action to an integrator to update display."""
            if hasattr(self, "Action"):
                raise ValueError("Action already added to an integrator.")
            self.Action = IntegrateObj.AddAction(Fn = self.Update, TimeFreq = TimeFreq
                                                 Name = "Visualizer")
    
        def DelAction(self, IntegrateObj):
            """Removes an action from an integrator."""
            IntegrateObj.DelAction(self.Action)
            del self.Action
    
