#!/usr/bin/env python


import re

import numpy as np
from var import *


MaxLen = 10000000



def Unique(l):
    return [x for (i,x) in enumerate(l) if not x in l[:i]]


def Uncomment(Source):
    l = []
    for line in Source.split("\n"):
        if "!" in line:
            line = line[:line.index("!")]
        l.append(line)
    return "\n".join(l)


def CombineLines(Source):
    """Combines any split lines into single ones."""
    p = re.compile(r"[\s]*\&[\s]*^[\s]*\&[\s]*", re.MULTILINE)
    if p.search(Source) is None:
        return Source
    else:
        return p.sub(" ", Source)
        
def CodeStrip(Source):
    """Strips code of beginning and ending whitespace, maintaining indents."""
    #convert tabs
    Source = Source.replace("\t", " "*4)
    #remove trailing space
    Source = Source.rstrip()
    #check if there is blank space at the beginning
    while "\n" in Source:
        FirstLine = Source.split("\n")[0]
        if len(FirstLine.strip()):
            return Source
        else:
            Source = Source[len(FirstLine)+1:]
    return Source

def SplitPoints(Line):
    """Returns points at which Line can be split."""
    #first aggregate operators and spaces
    l = []
    for (i, ch) in enumerate(Line):
        if ch in " +-,":
            l.append(i)
    #now find incides that can't be split
    p = re.compile(r"[\w][\w\d]*\(", re.IGNORECASE)
    for mobj in p.finditer(Line):
        Start = mobj.start()
        s = MatchClosingParen(Line[Start:])
        Stop = Start + len(s)
        #filter out inside of indices
        l = [i for i in l if i < Start or i >= Stop]
        #add commas back in
        for (k,ch) in enumerate(s):
            if ch == ",":
                l.append(Start + k + 1)
    l.sort()
    return l

def SplitLine(Line, SplitLen = MaxLen, Indent = 0, Contin = False):
    """Splits a line into a list of multiple lines."""
    Line = Line.rstrip()
    #check the length
    if len(Line) + Indent + Contin*4 < SplitLen:
        return [" "*Indent + "  & "*Contin + Line]
    #find the first non-space character
    l = 0
    for Char in Line:
        if Char == " ":
            l += 1
        else:
            break
    #find the space at which to split
    Points = SplitPoints(Line)
    Points.reverse()
    r = 0
    for Point in Points:
        if Point + Indent + Contin*4 < SplitLen - 2:
            r = Point
            break
    if r == 0 or r <= l:
        #don't use indent
        if len(Line) + Contin*4 < SplitLen:
            return ["  & "*Contin + Line]
        r = 0
        Indent = 0
        for Point in Points:
            if Point + Contin*4 < SplitLen - 2:
                r = Point
                break
        if r == 0:
            raise StandardError("Could not find a place to split line:\n%s" % Line)
    #add the first line
    Lines = [" "*Indent + "  & "*Contin + Line[:r] + " &"]
    #find subsequent lines by recursion
    Line = Line[r:].strip()
    if Indent == 0 and not Contin:
        Indent = l
    NewLines = SplitLine(Line, SplitLen, Indent = Indent, Contin = True)
    Lines.extend(NewLines)
    return Lines
    

def I(Source, Indent = 0, MaxLen = MaxLen):
    """Normalizes Fortran code."""
    if Source.strip() == "": return "\n"
    Source = CombineLines(Source).strip("\n")
    Lines = Source.split("\n")
    Lines = [" "*Indent + Line for Line in Lines]
    if MaxLen is None:
        return "\n".join(Lines) + "\n"
    else:       
        NewLines = []
        for Line in Lines:
            NewLines.extend(SplitLine(Line, MaxLen))
        return "\n".join(NewLines) + "\n"


def AllTokens(Source):
    """Returns a list of all tokens in code."""
    p = re.compile(r"(?<!\w)\w*[a-zA-Z]+\w*(?!\w)",
                   re.IGNORECASE)
    return Unique([m.group(0) for m in p.finditer(Uncomment(Source))])
        

def HasToken(Source, Token):
    """True if code has token."""
    p = re.compile(r"(?<!\w)"
                   + re.escape(Token)
                   + r"(?!\w)",
                   re.IGNORECASE)
    return (not p.search(Source) is None)

def ReplaceToken(Source, Token, Val):
    """Replaces Token in Source with Val."""
    p = re.compile(r"(?<!\w)"
                   + re.escape(Token)
                   + r"(?!\w)",
                   re.IGNORECASE)
    return p.sub(Val, Source)

def HasAssign(Source, VarL, VarR = None):
    """True if code makes an assignment to VarL."""
    if VarR is None:
        p = re.compile(r"(?<!\w)"
                       + re.escape(VarL)
                       + r"(?!\w).*(?<!\=)\=(?!\=)",
                       re.IGNORECASE)
        return (not p.search(Source) is None)
    else:
        p = re.compile(r"(?<!\w)"
                       + re.escape(VarL)
                       + r"(?!\w).*(?<!\=)\=(?!\=).*(?<!\w)"
                       + re.escape(VarR)
                       + r"(?!\w)",
                       re.IGNORECASE)
        return (not p.search(Source) is None)


def GetAssign(Source, Var):
    """Returns the expression to which Var is assigned, or zero length string if None."""
    p = re.compile(r"^[\s]*(?<!\w)("
                   + re.escape(Var)
                   + r"(?!\w).*)(?<!\=)\=(?!\=)(.*)$",
                   re.IGNORECASE | re.MULTILINE)
    ret = ""
    for m in p.finditer(Source):
        ret = m.group(2).strip()
    return ret


def ExtractAssign(Source, Var):
    """Returns the expression to which Var is assigned, or zero length string if None,
and removes from Source."""
    NewSource = ""
    Assign = ""
    for line in Source.split("\n"):
        s = GetAssign(line, Var)
        if len(s):
            Assign = s
        else:
            NewSource += line + "\n"
    return NewSource, Assign


def AugmentAssign(Source, Var, Coef = None):
    """Changes a = b to a = a + b."""
    p = re.compile(r"(?<!\w)("
                   + re.escape(Var)
                   + r"(?!\w).*)(?<!\=)\=(?!\=)(.*)",
                   re.IGNORECASE)
    def fn(mobj):
        s1 = mobj.group(1).strip()
        s2 = mobj.group(2).strip()
        if s2.startswith("-"):
            s2 = s2[1:]
            if Coef is None:
                return "%s = %s - %s" % (s1, s1, s2)
            else:
                return "%s = %s + (%s) * (-%s)" % (s1, s1, Coef, s2)
        else:
            if s2.startswith("+"): s2 = s2[1:]
            if Coef is None:
                return "%s = %s + %s" % (s1, s1, s2)
            else:
                return "%s = %s + (%s) * (%s)" % (s1, s1, Coef, s2)
    return p.sub(fn, Source)


def MatchParen(Source):
    """Returns code subset with a balanced number of parenthesis."""
    if len(Source) == 0: return Source
    ind = -1
    Diff = 0
    for (i,x) in enumerate(Source):
        if x == '(':
            Diff += 1
        elif x == ')':
            Diff -= 1
        if Diff < 0:
            break
        elif Diff == 0:
            ind = i
    return Source[:ind+1]

def MatchClosingParen(Source):
    """Returns code trimmed right to have balanced number of parenthesis."""
    if len(Source) == 0: return Source
    Diff = 0
    Found = False
    ind = -1
    for (i,x) in enumerate(Source):
        if x == '(':
            Diff += 1
            Found = True
        elif x == ')':
            Diff -= 1
        if Diff == 0 and Found:
            ind = i
            break
    return Source[:ind+1]
    
def ExtractParen(Source):
    """Returns anything inside the largest parenthetical enclosure."""
    if len(Source) == 0: return ""
    lind = Source.find('(')
    if lind < 0: return ""
    rind = Source.find(')')
    if rind < 0 or rind < lind: return ""
    return Source[rind:lind]
    
def SplitParen(Source):
    """Takes 'a(b)' and returns a tuple 'a', 'b'"""
    if len(Source) == 0: return "", ""
    lind = Source.find('(')
    rind = Source.find(')')
    if rind < 0 and lind < 0:
        return Source, ""
    elif rind > 0 and lind > 0 and rind > lind:
        return Source[:lind], Source[lind+1:rind]
    else:
        raise ValueError("Could not split out parenthetical enclosure in:\n%s" % Source)
    return Source[rind:lind]


def TrimParen(Source):
    """Removes opening and closing characters."""
    Source = Source.strip()
    if len(Source) == 0: return ""
    while Source[0] == '(' and Source[-1] == ')':
        Source = Source[1:-1]
        Source = Source.strip()
    return Source

def AddParen(Source):
    """Adds parenthesis if needed to Source."""
    m = re.match("[\w]*", Source)
    if m is None:
        return "(%s)" % Source
    else:
        if len(m.group(0)) == len(Source):
            return Source
        else:
            return "(%s)" % Source    


def GetInd(Source):
    """Returns a list of string indices."""
    Source = MatchParen(Source)
    Source = TrimParen(Source)
    return [x.strip() for x in Source.split(',')]


def Eval(Source):
    """Evaluates any integer additions / subtractions / multiplications if possible."""
    Source = Source.strip()
    if len(Source) == 0: return ""
    p1 = re.compile(r"([\s\+\-\(]|^)([\d]+[\s]*[\+\-\*][\s]*[\d]+)([\s\+\-\)]|$)")
    p2 = re.compile(r"\([\s]*([\d]+)[\s]*\)")
    def fn1(mobj):
        s = mobj.group(0)
        t = mobj.group(1) + str(eval(mobj.group(2))) + mobj.group(3)
        return t
    def fn2(mobj):
        return mobj.group(1)
    while not p1.search(Source) is None:
        Source = p1.sub(fn1, Source)
        Source = p2.sub(fn2, Source)
    return Source



def ReplaceVarInd(Source, ReplaceFn):
    NewSource = Source
    p = re.compile(r"(?<![a-zA-Z0-9_])"
                   + re.escape(OldVar)
                   + r"\(",
                   re.IGNORECASE)
    for m in p.finditer(Source):
        IndSpec = MatchParen(Source[m.end(0):])
        #check for only one index
        Indices = GetInd(IndSpec)
        NewVar = ReplaceFn(Indices)
        #make the replacement
        Token = "%s(%s)" % (OldVar, IndSpec)
        NewSource = ReplaceToken(NewSource, Token, NewVar)
    return NewSource    
    

def ReplaceInd(Source, Strides = (1,), StartInd = 1, ShiftInd = None):
    if ShiftInd is None:
        ShiftInd = [0]*len(Strides)
    if not len(ShiftInd) == len(Strides):
        raise ValueError("ShiftInd and Strides must be same length.")
    Indices = GetInd(Source)
    if not len(Indices) == len(Strides):
        raise IndexError, "Found %d indices but found %d elements in %s" % \
              (len(Strides), len(Indices), Source)
    s = []
    for (Shift, Stride, ind) in zip(ShiftInd, Strides, Indices):
        if Stride == 1:
            s.append(Eval(ind))
            StartInd += Shift
        else:
            sInd = AddParen(Eval(ind))
            s.append(Eval("%d*%s" % (Stride, sInd)))
            StartInd += Stride * Shift
    if StartInd > 0:
        s = ["%s" % StartInd] + s
    s = " + ".join(s)
    s = TrimParen(Eval(s))
    return s
    

def ReplaceVarMaster(Source, OldVar, NewVar,
                     Strides = (1,), StartInd = 1, ShiftInd = None,
                     Len = None):
    """Replaces an indexed variable with a new name and index
in a master array."""
    NewSource = Source
    #first do parentheses
    p = re.compile(r"(?<![a-zA-Z0-9_])"
                   + re.escape(OldVar)
                   + r"\(",
                   re.IGNORECASE)
    for m in p.finditer(Source):
        IndSpec = MatchParen(Source[m.end(0):])
        #get index expression for 1D array
        IndExpr = ReplaceInd(IndSpec, Strides, StartInd, ShiftInd)
        #make the replacement
        Token = "%s(%s)" % (OldVar, IndSpec)
        Val = "%s(%s)" % (NewVar, IndExpr)
        NewSource = ReplaceToken(NewSource, Token, Val)
    #next do non-indexed
    p = re.compile(r"(?<![a-zA-Z0-9_])"
                   + re.escape(OldVar)
                   + r"(?![a-zA-Z0-9_\(])",
                   re.IGNORECASE)
    if ShiftInd is None and len(Strides) == 1:
        ShiftInd = [0]
    for m in p.finditer(Source):
        if len(Strides) > 1:
            raise ValueError("Must specify explicit array indices for variable %s" % OldVar)
        #make the replacement
        Token = "%s" % OldVar
        if Len == 1:
            Val = "%s(%d)" % (NewVar, StartInd + ShiftInd[0])
        else:
            Val = "%s(%d:%d)" % (NewVar, StartInd + ShiftInd[0], StartInd+ShiftInd[0]+Len-1)
        NewSource = ReplaceToken(NewSource, Token, Val)
    return NewSource

def GetMasterIndList(Source, Var,
                     Strides = (1,), StartInd = 1, ShiftInd = (0,)):
    """Returns all indices referenced to a Master array."""
    p = re.compile(r"(?<![a-zA-Z0-9_])"
                   + re.escape(Var)
                   + r"\(",
                   re.IGNORECASE)
    IndList = []
    for m in p.finditer(Source):
        IndSpec = MatchParen(Source[m.end(0):])
        #get index expression for 1D array
        IndExpr = ReplaceInd(IndSpec, Strides, StartInd, ShiftInd)
        #add to list
        if not IndExpr in IndList:
            IndList.append(IndExpr)
    return IndList    
        

def IfBracket(Source, Condition, Indent = 4):
    """Places Source within an if bracket."""
    if Condition:
        if Source.strip():
            s = I("if (%s) then" % Condition)
            s += I(Source, Indent)
            s += I("end if")
            return s
        else:
            return Source
    else:
        return Source
        
def BooleanAnd(Source1, Source2, Indent = 4):
    """Places Source within an if bracket."""
    if Source1:
        if Source2:
            return "(%s) .and. (%s)" % (Source1, Source2)
        else:
            return Source1
    else:
        if Source2:
            return Source2
        else:
            return ""
            
def BooleanOr(Source1, Source2, Indent = 4):
    """Places Source within an if bracket."""
    if Source1:
        if Source2:
            return "(%s) .or. (%s)" % (Source1, Source2)
        else:
            return Source1
    else:
        if Source2:
            return Source2
        else:
            return ""            
        
def npStrides(Strides, dtype):
    itemsize = np.dtype(dtype).itemsize
    return tuple([itemsize*x for x in Strides])

def Strides(a):
    return tuple([x/a.itemsize for x in a.strides])
    
 
    

        
    

        