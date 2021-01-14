#!/usr/bin/env python

import re


Order = 'F'

Types = {"float":"real(8)", "int":"integer", "bool":"logical", "complex":"complex*16"}


def VarSplit(Var, FortranType = False):
    """Splits an argument into a common block, type, name, shape, and value."""
    Var = Var.strip()
    rePTypes = "|".join(Types.keys())
    if Var.startswith("/"):
        #common block spec; this must be a def
        if "=" in Var:
            raise ValueError("Cannot specify a common block with a mapped variable:\n%s" % Var)
        i = Var.find("/",1)
        if i < 0:
            raise ValueError("Could not find closing '/' in common block for variable:\n%s" % Var)
        CommonBlock = Var[1:i].strip()
        Var = Var[i+1:]
        Val = ""
    elif "=" in Var:
        Var, Val = Var.split("=")
        Var = Var.strip()
        Val = Val.strip()
        CommonBlock = ""
    else:
        CommonBlock = ""
        Val = ""
    if not " " in Var and not "(" in Var:
        #this is an external without a type (or an error)
        return "", "", Var, "", ""
    p = re.compile(r"(?P<type>" + rePTypes + r")\s*"
                   + r"(?P<name>\w*)"
                   + r"(?P<shape>\(.*\))?\s*",
                   re.IGNORECASE | re.MULTILINE)
    m = p.search(Var)
    if m is None:
        raise TypeError("Could not parse variable spec %s" % Var)
    Type = m.group("type")
    Name = m.group("name")
    Shape = m.group("shape")
    if Shape is None: Shape = ""
    #process the types
    if FortranType:
        Type = Types[Type.lower()]
    #proces the shape
    if len(Shape):
        #get the dimensions
        Shape = Shape[1:-1]
        ParseShape = []
        for d in Shape.split(","):
            d = d.strip()
            if ":" in d:
                ParseShape.append(d)
            else:
                ParseShape.append("0:%s-1" % d)
        Shape = ", ".join(ParseShape)
    return CommonBlock, Type, Name, Shape, Val
    

def VarJoin(Type, Name, Shape, Val, CBlock = ""):
    """Joins variable specs into a single argument."""
    if len(Val):
        if len(CBlock):
            raise ValueError("Cannot specify a common block with a mapped variable:\n%s" % Name)
        return "%s %s(%s) = %s" % (Type, Name, Shape, Val)
    elif len(CBlock):
        return "%s: %s %s(%s)" % (CBlock, Type, Name, Shape)
    else:       
        return "%s %s(%s)" % (Type, Name, Shape)


def VarType(Var, FortranType = False):
    """Returns the type of an argument."""
    CBlock, Type, Name, Shape, Val = VarSplit(Var, FortranType = FortranType)
    return Type


def VarName(Var):
    """Returns the type of an argument."""
    CBlock, Type, Name, Shape, Val = VarSplit(Var)
    return Name    


def VarShape(Var):
    """Returns variables used in a dimension specification."""
    CBlock, Type, Name, Shape, Val = VarSplit(Var)
    if not Shape: return []
    p = re.compile(r"(?<!\w)\w*[a-zA-Z]+\w*(?!\w)",
                   re.IGNORECASE)
    ShapeVars = []
    for m in p.finditer(Shape):
        if not m.group(0) in ShapeVars:
            ShapeVars.append(m.group(0))
    ShapeVars = ["int " + d.strip() for d in ShapeVars]
    return ShapeVars


def VarVal(Var):
    """Returns the value of an argument."""
    CBlock, Type, Name, Shape, Val = VarSplit(Var)
    return Val    

        
def VarNorm(Var):
    """Returns a normalized version with dimension ranges explicit."""
    CBlock, Type, Name, Shape, Val = VarSplit(Var)
    if Shape:
        s = "%s %s(%s)" % (Type, Name, Shape)
    else:
        s = "%s %s" % (Type, Name)
    if CBlock:
        if Val:
            raise ValueError("Cannot specify a common block with a mapped variable:\n%s" % Name)
        s = "/%s/" % CBlock + s        
    elif Val:
        s += " = %s" % Val
    return s
    
def IsAllocatable(Shape):
    """Returns True if variable is defined like Var(:) or Var(:,:), etc"""
    Shape = Shape.replace(" ","")
    return Shape == ":" or ":," in Shape 
    
def IsArray(Var):
    """Returns True if Var has a shape."""
    CBlock, Type, Name, Shape, Val = VarSplit(Var)
    return len(Shape) > 0


def GetVarType(Obj):
    """Returns a string type for Obj and IsArray""" 
    import numpy as np
    IsArray = False
    if isinstance(Obj, float):
        Type = "float"
    elif isinstance(Obj, complex):
        Type = "complex"
    elif isinstance(Obj, bool):
        Type = "bool"
    elif isinstance(Obj, int):
        #this must come after bool testing since isinstance(True, int) = True
        Type = "int"        
    elif type(Obj) is np.ndarray:
        IsArray = True
        if Obj.dtype is np.dtype(float):
            Type = "float"
        elif Obj.dtype is np.dtype(complex):
            Type = "complex"
        elif Obj.dtype is np.dtype(int):
            Type = "int"
        elif Obj.dtype is np.dtype(bool):
            Type = "bool"            
        else:
            raise TypeError("Do not recognize type of %s." % VarName)
    else:
        raise TypeError("Do not recognize type of %s." % VarName)
    return Type, IsArray


def TypeObject(Type):
    """Returns a function to convert to Python types."""
    if Type == "int":
        return int
    elif Type == "float":
        return float
    elif Type == "complex":
        return complex    
    elif Type == "bool":
        return bool
    else:
        raise ValueError("Unsupported type.")

        
def GetVarSpec(VarName, Obj, IsConstant = False, ShapeStr = None):
    """Makes a variable specification from an actual Python object.""" 
    from code import ReplaceToken
    Type, IsArray = GetVarType(Obj)
    s = "%s %s" % (Type, VarName)
    if IsArray:
        if ShapeStr:
            s += "(" + ShapeStr + ")"
        else:
            s += "(" + ", ".join(["0:%d-1" % x for x in Obj.shape]) + ")"
        if IsConstant:
            s += " = "
            l = len(Obj.shape)
            if l == 1:
                s += "(/ " + ", ".join([str(x) for x in Obj]) + " /)"
            elif l == 2:
                #have to do some mangling here because fortran reshape differs from numpy one
                a = Obj.T.flatten()
                sdata = ",".join([str(x) for x in a])
                s += "reshape( (/ %s /) , (/%d, %d/) )" % (sdata, Obj.shape[0], Obj.shape[1])
            else:
                raise ValueError("I do not know how to process an array with >= 3dimensions (not yet implemented).")
    elif IsConstant:
        s += " = " + str(Obj)
    s = ReplaceToken(s, "True", "1")
    s = ReplaceToken(s, "False", "0")
    return s
    

