#!/usr/bin/env python

### Simple trajectory format.
### coded by MSS

import base


class Mapped(base.TrajClass):
    """Wrapper for a mapping function applied to a trajectory."""

    def __init__(self, Traj, PosMap, Sys = None, Mask = None,
                  AtomNames = None, BoxL = None):
        """Wrapper for a mapping function applied to a trajectory."""
        import copy
        base.TrajClass.__init__(self, Sys=Sys, Mask=Mask)
        self.__PosMap = PosMap
        self.__Traj = Traj
        if AtomNames is None:
            self.AtomNames = [Map.Atom2Name for Map in PosMap]
        else:
            self.AtomNames = AtomNames
        self.BoxL = copy.copy(BoxL)
        self.NAtom = len(self.AtomNames)
        self.NFrameAbs = self.__Traj.NFrame
        self.SetSlice()

    def Close(self):
        self.__Traj.Close()

    def ParseFrame(self, IndexAbs):
        Pos1 = self.__Traj.Get(IndexAbs)        
        self.FrameData = self.__Traj.FrameData.copy()
        if "BoxL" in self.FrameData:
            BoxL = self.FrameData["BoxL"]
        elif not self.BoxL is None:
            BoxL = self.BoxL
        elif not self.Sys is None:
            BoxL = self.Sys.BoxL
        else:
            raise ValueError("Cannot find a specification for BoxL.")
        Pos2 = self.__PosMap.Map(Pos1, BoxL)
        return Pos2

    
        