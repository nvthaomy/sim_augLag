#/usr/bin/env python

### Class definitions for parameter arrays in SIM suite.
### coded by MSS

import numpy as np

import sim.fortran as fortran



def npStrides(Strides, dtype):
    itemsize = np.dtype(dtype).itemsize
    return tuple([itemsize*x for x in Strides])

def Strides(a):
    return tuple([x/a.itemsize for x in a.strides])

def LinkArray(Array, Master, StartInd):
    n = len(Array)
    Master[StartInd:StartInd + n] = Array
    Array.data = Master[StartInd:StartInd + n].data

def MakeMasterArray(ArrayList):
    LinkedArrays = []
    n = sum([len(a) for a in ArrayList])
    dtype = ArrayList[0].dtype
    Master = np.zeros(n, dtype=dtype)
    Ind = 0
    for a in ArrayList:
        LinkArray(a, Master, Ind)
        LinkedArrays.append((a, Ind))
        Ind += len(a)
    return Master, LinkedArrays

   

def MasterSizes(PGlobal, PList):
    """Returns the sizes of the master DUParam and DDUParam arrays."""
    NList = np.sum([P.Param.DDU.N for P in PList]) 
    N2List = np.sum([P.Param.DDU.N2 for P in PList])
    NGlobal = PGlobal.DDU.N
    N2Global = PGlobal.DDU.N2
    N2Cross = NList * NGlobal
    N = NList + NGlobal
    N2 = N2Global + N2List + N2Cross
    return N, N2


class Hessian(object):
    """Generic class for Hessian matrices;
stores data in a 1D array."""

    def __init__(self, N):
        self.Array = np.zeros(shape=N*N, dtype=float)
        self.N = N
        self.N2 = N*N
        Matrix = self.Array.reshape((self.N, self.N),
                                    order = fortran.Order)
        self.Strides = Strides(Matrix)

    def __getattr__(self, name):
        if name == "Matrix":
            Matrix = self.Array.reshape((self.N, self.N),
                                        order = fortran.Order)
            return Matrix + Matrix.T - np.diag(Matrix.diagonal())
        else:
            raise AttributeError(name)

    def __setitem__(self, key, val):
        if not type(key) is tuple:
            raise IndexError("Key must be integer tuple.")
        i,j = key
        k = i + j * self.N
        self.Array[k] = val

    def __getitem__(self, key):
        if not type(key) is tuple:
            raise IndexError("Key must be integer tuple.")
        i,j = key
        k1 = i + j * self.N
        k2 = j + i * self.N
        return self.Array[k1] + self.Array[k2]           


class CrossHessian(object):
    """Generic class for Hessian matrices between two sets of vars;
stores data in a 1D array."""

    def __init__(self, Ni, Nj):
        self.Array = np.zeros(shape=Ni*Nj, dtype=float)
        self.Ni = Ni
        self.Nj = Nj
        self.N2 = Ni*Nj
        self.Matrix = self.Array.reshape((self.Ni, self.Nj),
                                         order = fortran.Order)
        self.Strides = Strides(self.Matrix)

    def __setitem__(self, key, val):
        if not type(key) is tuple:
            raise IndexError("Key must be integer tuple.")
        self.Matrix[key] = val

    def __getitem__(self, key):
        if not type(key) is tuple:
            raise IndexError("Key must be integer tuple.")
        return self.Matrix[key]         


class MasterHessian(object):
    """Generic class for master Hessian matrices;
stores data in a 1D array."""

    ModuleVars = ["Array"]

    def __init__(self, HessianGlobal, HessianList, dtype = float):
        self.NList = np.sum([a.N for a in HessianList]) 
        self.N2List = np.sum([a.N2 for a in HessianList])
        self.NGlobal = HessianGlobal.N
        self.N2Global = HessianGlobal.N2
        self.N2Cross = self.NList * self.NGlobal
        self.N = self.NGlobal + self.NList 
        self.N2 = self.N2Global + self.N2List + self.N2Cross
        self.Array = np.zeros(shape=self.N2, dtype=float)
        self.__Range = []
        self.Strides = []
        self.LinkedArrays = []
        Ind1, Ind2 = 0, 0
        #make the sub-hessians
        for h in [HessianGlobal] + HessianList:
            LinkArray(h.Array, self.Array, Ind2)
            self.LinkedArrays.append((h.Array, Ind2))
            self.Strides += [h.Strides]
            self.__Range.append((Ind2, Ind1, Ind1 + h.N, Ind1, Ind1 + h.N))
            Ind1 += h.N
            Ind2 += h.N2
        #make the cross-hessian between globals and others
        Ind1 = self.NGlobal
        for h in HessianList:
            Cross = CrossHessian(HessianGlobal.N, h.N)
            LinkArray(Cross.Array, self.Array, Ind2)
            self.LinkedArrays.append((Cross.Array, Ind2))
            self.Strides += [Cross.Strides]
            self.__Range.append((Ind2, 0, self.NGlobal, Ind1, Ind1 + h.N))
            Ind1 += h.N
            Ind2 += Cross.N2
            
    def RelinkArrays(self):
        """Relinks arrays to the master array."""
        for (Array, Ind) in self.LinkedArrays:
            LinkArray(Array, self.Array, Ind)
                 
    def GetInd(self, i, j):
        """Returns the 1D index for a 2D index."""
        for (k, (Ind2, iStart, iStop, jStart, jStop)) in enumerate(self.__Range):
            if i >= iStart and i < iStop and j >= iStart and j < jStop:
                Strides = self.Strides[k]
                return Ind2 + Strides[0]*(i - iStart) + Strides[1]*(j - jStart)
        raise ValueError("i,j not from same sub-matrix.")

    def MakeMatrix(self, Array = None):
        """Returns a matrix version of data."""
        if Array is None:
            Array = self.Array
        else:
            if not Array.shape == self.Array.shape:
                raise TypeError("Shape of input array does not match internal one.")
        ret = np.zeros((self.N,self.N), dtype=Array.dtype)
        Start2 = 0
        for (k, (Ind2, iStart, iStop, jStart, jStop)) in enumerate(self.__Range):
            Ni = iStop - iStart
            Nj = jStop - jStart
            strides = npStrides(self.Strides[k], Array.dtype)
            subm = np.ndarray(shape=(Ni, Nj), dtype=Array.dtype,
                              buffer = Array.data,
                              offset = Start2 * Array.itemsize,
                              strides = strides)
            ret[iStart:iStop, jStart:jStop] = subm
            Start2 += Ni * Nj
        ret = ret + ret.T - np.diag(ret.diagonal())
        return ret    

    def __getattr__(self, name):
        if name == "Matrix":
            return self.MakeMatrix()
        else:
            raise AttributeError(name)

     


class Parray(object):
    """Generic class for describing and containing parameters."""
    Array1Attr = ["Fixed", "Min", "Max", "SoftMin", "SoftMax", "DU", "DW",
                  "MaxChange", "Scale"]
    Array1Vars = Array1Attr + ["Val"] 
    Array2Vars = ["DDU", "DDW"]
    ArrayVars = Array1Vars + Array2Vars

    def __init__(self, Name = ""):
        self.Maps = {}
        self.Name = Name
        self.Names = []
        self.N = 0
        self.N2 = 0
        self.Ind1 = None
        self.Ind2 = None
        self.Ind2Global = None
        self.Fixed = np.zeros(0, dtype=bool)
        self.Min = np.ones(0, dtype=float) * -1.e300
        self.Max = np.ones(0, dtype=float) * 1.e300
        self.SoftMin = np.ones(0, dtype=float) * -1.e300
        self.SoftMax = np.ones(0, dtype=float) * 1.e300
        for v in Parray.Array1Vars[5:]:
            self.__dict__[v] = np.zeros(0, dtype=float)
        for v in Parray.Array2Vars:
            self.__dict__[v] = Hessian(0)
        self.Const = np.zeros(0, dtype=float)
      
    def Add(self, Name, N = 1,
            Value = None, Fixed = None,
            Min = None, Max = None,
            SoftMin = None, SoftMax = None,
            MaxChange = 1.e300, Scale = 1.):
        if Name in self.__dict__:
            raise ValueError("Cannot set map with protected name %s" % Name)
        Start, Stop = self.N, self.N + N
        self.N += N
        self.N2 = self.N*self.N
        self.Maps[Name] = (Start, Stop)
        if N > 1:
            for i in range(N):
                self.Names.append("%s(%d)" % (Name, i))
        elif N == 1:
            self.Names.append("%s" % Name)
        #extend arrays
        for v in self.Array1Vars:
            var = getattr(self, v)
            delattr(self, v)
            setattr(self, v, np.concatenate((var, np.zeros(N,var.dtype))))
        for v in self.Array2Vars:
            var = getattr(self, v)
            delattr(self, v)
            setattr(self, v, Hessian(self.N))
        #place values
        if not Value is None:
            self.Val[Start:Stop] = Value
        if not Fixed is None:
            self.Fixed[Start:Stop] = Fixed
        self.Min[Start:Stop] = -1.e300
        self.SoftMin[Start:Stop] = -1.e300
        self.Max[Start:Stop] = 1.e300
        self.SoftMax[Start:Stop] = 1.e300
        if not Min is None:
            self.Min[Start:Stop] = Min
            self.SoftMin[Start:Stop] = Min
        if not Max is None:
            self.Max[Start:Stop] = Max
            self.SoftMax[Start:Stop] = Max
        if not SoftMin is None:
            self.SoftMin[Start:Stop] = SoftMin
        if not SoftMax is None:
            self.SoftMax[Start:Stop] = SoftMax
        self.MaxChange[Start:Stop] = MaxChange
        self.Scale[Start:Stop] = Scale

    def GetInd(self, name):
        if name in self.Maps:
            return self.Maps[name][0]
        else:
            raise AttributeError(name)    

    def __getattr__(self, name):
        if name in self.Maps:
            Start, Stop = self.Maps[name]
            N = Stop - Start
            return ParrayMap(shape = N, dtype = float,
                             buffer = self.Val.data,
                             offset = Start * self.Val.itemsize,
                             Parray = self, Start = Start,
                             Ind1 = self.Ind1)           
        else:
            raise AttributeError(name)

    def __setattr__(self, name, val):
        if name in self.ArrayVars and name in self.__dict__:
            getattr(self, name)[:] = val
        elif name in self.__dict__.get("Maps", []):
            i,j = self.Maps[name]
            self.Val[i:j] = val
        else:
            object.__setattr__(self, name, val)
           

class ParrayMap(np.ndarray):
    """Class for containing mapped parameters."""
    def __new__(self, shape, dtype=float, buffer=None, offset=0,
                strides=None, order=None,
                Parray=None, Start=None, Ind1=None):
        obj = np.ndarray.__new__(self, shape=shape, dtype=dtype,
                                 buffer=buffer, offset=offset,
                                 strides=strides, order=order)
        n = len(obj)
        obj.Fixed = np.zeros(n, dtype=bool)
        obj.Min = np.ones(n, dtype=float) * -1.e300
        obj.Max = np.ones(n, dtype=float) * 1.e300
        obj.SoftMin = np.ones(n, dtype=float) * -1.e300
        obj.SoftMax = np.ones(n, dtype=float) * 1.e300
        if Ind1 is None:
            obj.Ind1 = None
        else:
            obj.Ind1 = Ind1 + Start
        for v in Parray.Array1Attr[5:]:
            setattr(obj, v, np.zeros(n, dtype=float))
        if not Parray is None and not Start is None:
            for v in Parray.Array1Attr:
                pvar = getattr(Parray, v)[Start:Start+n]
                var = getattr(obj, v)
                var.data = pvar.data
        return obj

    def __setattr__(self, name, val):
        if name in Parray.Array1Attr and hasattr(self, name):
            getattr(self, name)[:] = val
        else:
            np.ndarray.__setattr__(self, name, val)

           
class MasterParray(object):
    """Generic master class for describing and containing parameters."""
    
    
    def __init__(self, ParrayGlobal, ParrayList, Name = ""):
        self.Name = Name
        self.Names = []
        self.N = 0
        self.N2 = 0
        self.LinkedArrays = {}
        AllParrays = [ParrayGlobal] + ParrayList
        for a in AllParrays:
            self.Names += ["%s:%s" % (a.Name, nm) for nm in a.Names]
            a.Ind1 = self.N
            a.Ind2 = self.N2
            self.N += a.N
            self.N2 += a.N2
        for a in ParrayList:
            a.Ind2Global = self.N2
            self.N2 += ParrayGlobal.N * a.N
        for Var in Parray.Array1Vars:
            ArrayList = [getattr(a, Var) for a in AllParrays]
            Master, LinkedArrays = MakeMasterArray(ArrayList)
            self.LinkedArrays[Var] = LinkedArrays
            setattr(self, Var, Master)
        for Var in Parray.Array2Vars:
            HessianGlobal = getattr(ParrayGlobal, Var)
            HessianList = [getattr(a, Var) for a in ParrayList]
            Master = MasterHessian(HessianGlobal, HessianList)
            setattr(self, Var, Master)
            #update strides
            ParrayGlobal.Strides2 = Master.Strides[0]
            for (i,a) in enumerate(ParrayList):
                a.Strides2 = Master.Strides[i+1]
                a.Strides2Global = Master.Strides[i+1+len(ParrayList)]

    def RelinkArrays(self):
        """Relinks all sub arrays to masters."""
        for Var in Parray.Array1Vars:
            Master = getattr(self, Var)
            for (Array, Ind) in self.LinkedArrays[Var]:
                LinkArray(Array, Master, Ind)
        for Var in Parray.Array2Vars:
            Master = getattr(self, Var)
            Master.RelinkArrays()
    
    def __setattr__(self, name, val):
        if name in Parray.ArrayVars and name in self.__dict__:
            getattr(self, name)[:] = val
        else:
            object.__setattr__(self, name, val)
            
    def PreLoad(self, Lib):
        """Run before compiling."""
        self.NDParam = len(self.Val)
        self.NDDParam = len(self.DDU.Array)
        Lib.VarPreLoad(self, Vars = [("Val", "Param"), ("DU", "DUParam"), ("DW", "DWParam")],
                       Consts = ["NDParam", "NDDParam"])
        Lib.VarPreLoad(self.DDU, Vars = [("Array", "DDUParam")])
        Lib.VarPreLoad(self.DDW, Vars = [("Array", "DDWParam")])
        
    def PostLoad(self, Lib):
        """Run before compiling."""
        Lib.VarPostLoad(self)
        Lib.VarPostLoad(self.DDU)
        Lib.VarPostLoad(self.DDW)
        #relink arrays
        self.RelinkArrays()


            


    