#/usr/bin/env python

### Fortran code block generation

import copy
from var import VarName, VarShape, VarSplit, VarJoin
from code import HasToken, AllTokens, I, ReplaceToken, IfBracket, CodeStrip, SplitParen


class CodeBlock(object):
    
    def __init__(self, Source = None):
        """Initializes a new CodeBlock block."""
        self.Clear()
        if not Source is None:
            self.Add(Source)
        
    def copy(self):
        """Returns a copy of self."""
        return copy.deepcopy(self)
        
    def Copy(self):
        """Returns a copy of self."""
        return self.copy()
        
    def RawSource(self):
        """Returns the raw source for this object."""
        if not len(self.Lines):
            return ""
        s = []
        MaxLen = max([len(l) + i for (l,d,i) in self.Lines])
        for (l,d,i) in self.Lines:
            ThisLen = len(l) + i
            ThisLine = " "*i + l
            if type(d) is str:
                ThisLine += " "*(MaxLen - ThisLen) 
                ThisLine += "  || %s" % d
            elif type(d) is list:
                ThisLine += " "*(MaxLen - ThisLen)
                ThisLine += "  || %s" % (",".join(d))
            s.append(ThisLine) 
        return "\n".join(s)
        
    def Clear(self):
        """Clears contents of self."""
        self.Lines = []
        self.ProcessedSource = None
        self.Flags = set()
        
    def AddHolder(self, HolderName, Dep = None, Indent = 0):
        """Adds a holder token."""
        self.Add("[%s]" % HolderName, Dep, Indent)
        
    def HasHolder(self, HolderName):
        """Returns True if has a holder."""
        n = "[%s]" % HolderName
        m = "[%s(" % HolderName
        for (s,d,i) in self.Lines:
            if s == n:
                return True
            elif s.startswith(m):
                return True
        return False
            
    def Add(self, Source = "", Dep = None, Indent = 0):
        """Adds additional source."""
        if not type(Dep) in [type(None), str, list]:
            raise TypeError("Dep must be None, string, or list; found %s:\n%s" % (type(Dep), str(Dep)))
        Source = CodeStrip(Source)
        if len(Source):
            if type(Dep) is str:
                Dep = Dep.strip()
            elif type(Dep) is list:
                Dep = [x.strip() for x in Dep]
            UseLastDep = False
            LastDep = Dep
            for line in Source.split("\n"):
                #get indent
                lIndent = len(line) - len(line.lstrip())
                #look for specified dependencies
                if "||" in line:
                    #see if we are to continue using this dependency for subsequent lines
                    if "||+" in line:
                        line, lDep = line.split("||+")[:2]
                        UseLastDep = True
                    else:
                        line, lDep = line.split("||")[:2]
                        UseLastDep = False
                    #make a list of dependencies
                    lDep = lDep.split(",")
                    lDep = [x.strip() for x in lDep]
                    #add global dependencies
                    if type(Dep) is list:
                        lDep = lDep + Dep                      
                    elif type(Dep) is str:
                        lDep = lDep + [Dep]
                    #filter out duplicates
                    lDep = [x for (i,x) in enumerate(lDep) if not x in lDep[i+1:]]
                    #simplify
                    if len(lDep) == 0: 
                        lDep = None
                    elif len(lDep) == 1:
                        lDep = lDep[0]
                    #record last dependency
                    LastDep = lDep
                elif UseLastDep:
                    #see if we are to use dependencies
                    lDep = LastDep
                else:
                    lDep = Dep
                line = line.strip()
                self.Lines.append((line, lDep, Indent + lIndent))
        else:
            self.Break()
            
    def Break(self):
        """Adds a line feed to source code."""
        self.Lines.append(("", None, 0))
            
    def AddFlag(self, Flag):
        """Adds a dependency flag to be checked."""
        self.Flags.add(Flag)
            
    def AddCodeBlock(self, other):
        """Adds other to current source code."""
        if not isinstance(other, CodeBlock):
            raise TypeError("Cannot add; other is not of type CodeBlock")
        o = other.copy()
        self.Lines.extend(o.Lines)
        self.Flags.update(o.Flags)
        
    def ReplaceTokenSource(self, Token, Val):
        """Replaces Token with Val in Source lines."""
        self.Lines = [(ReplaceToken(s, Token, Val), d, i) for (s, d, i) in self.Lines]
                       
    def Indent(self, Indent):
        """Indents each line of source by bulk amount."""
        self.Lines = [(s, d, i+Indent) for (s,d,i) in self.Lines]
        self.Indent = 0

    def Process(self, HolderDict = {}):
        """Returns processed code based on dependencies."""
        #make a copy of flags
        Flags = self.Flags.copy()
        #add holder dictionary names to flags
        Flags.update(["[%s]" % k for (k,v) in HolderDict.items() if len(v.strip())])
        #initial source
        LastSource = ""
        #loop over and add dependencies as necessary
        while True:
            LastLen = 0
            NewSource = ""
            for (s, Dep, Indent) in self.Lines:
                if Dep is None:
                    DoAdd = True
                elif type(Dep) is str: 
                    DoAdd = HasToken(LastSource, Dep) or Dep in Flags
                elif type(Dep) is list:
                    DoAdd = any([HasToken(LastSource, d) or d in Flags for d in Dep])
                else:
                    raise TypeError("Dep must be None, string, or list.")
                if DoAdd:
                    if s.startswith("[") and s.endswith("]"):
                        #this is a holder block
                        #split out [blockname(var1, var2, ...)] into "blockname" and ["var1", "var2", ...]
                        key, replacements = SplitParen(s[1:-1])
                        if key in HolderDict:
                            s = HolderDict[key]
                            #replace assignments using "Var1 = ThisVar" to sub "ThisVar" for "Var" 
                            if len(replacements):
                                for l in replacements.split(","):
                                    OldToken, NewToken = l.split("=", 1)
                                    OldToken = OldToken.strip()
                                    NewToken = NewToken.strip()
                                    s = ReplaceToken(s, OldToken, NewToken)
                            s = CodeStrip(s)
                            ThisLen = len(s.strip())
                            if ThisLen:
                                NewSource = NewSource + I(s, Indent)
                                LastLen = ThisLen
                        else:
                            #don't truncate the line length because this is a holder name we've not yet replaced
                            s = CodeStrip(s)
                            NewSource = NewSource + I(s, Indent, MaxLen = 10000000)
                            LastLen = len(s)
                    else:
                        ThisLen = len(s.strip())
                        if ThisLen > 0 or LastLen > 0:
                            NewSource = NewSource + I(s, Indent)
                        LastLen = ThisLen
            NewSource = CodeStrip(NewSource)
            if len(NewSource) == len(LastSource):
                self.ProcessedSource = NewSource
                return NewSource
            else:
                LastSource = NewSource



class FortCode(object):
    
    def __init__(self, Source = "", Blocks = [], VarIndex = None, Comment = ""):
        """Initializes a FortCode object."""
        #Source is a list of code blocks with dependency lists
        self.ProcessedSource = None
        self.InputMaps = []
        self.OutputMaps = []
        self.Defs = []
        self.Externals = []
        self.Reserved = []
        self.Blocks = {}
        self.CurBlock = "Main"
        self.VarIndex = VarIndex
        for x in ["Main"] + Blocks:
            self.NewBlock(x)
        self.Add(Source)
        self.Comment = Comment
        self.IfCondition = None

    def HasDef(self, Var):
        """Returns True if variable name in Var is already defined."""
        nVar = VarName(Var).lower()
        for Def in self.Defs:
            if nVar == VarName(Def).lower():
                return True
        return False
        
    def HasBlock(self, Block):
        """Returns True if Block is present."""
        return Block in self.Blocks    
        
    def ReferencesBlock(self, BlockName):
        """Returns True if Block is present."""
        return any([BlockName in Block.RawSource() for Block in self.Blocks.values()])              
            
    def NewUniqueDef(self, Var, ForceSuffix = False):
        CBlock, Type, Name, Shape, Val = VarSplit(Var)
        if self.HasDef(Var) or ForceSuffix:
            i = 0
            while True:
                NewName = "%s_%d" % (Name, i)
                NewVar = VarJoin(Type, NewName, Shape, Val, CBlock)
                if not self.HasDef(Var):
                    self.Defs.append(Var)
                    return NewVar
                i += 1
        else:
            self.Defs.append(Var)
            return Var
        
    def NewBlock(self, Name):
        """Adds a new source block."""
        if Name in self.Blocks:
            raise ValueError("Already have object %s" % Name)
        cb = CodeBlock()
        self.Blocks[Name] = cb
        return cb
    
    def DelBlock(self, Name):
        """Deletes a block."""
        if Name in self.Blocks:
            del self.Blocks[Name]
        else:
            raise KeyError(Name)
        
    def SetCurrent(self, Block):
        """Sets the default source block for adding."""
        if not Block in self.Blocks:
            raise ValueError("%s not found in Blocks." % Block)
        self.CurBlock = Block
        
    def copy(self):
        """Returns a copy of self."""
        return copy.deepcopy(self)
        
    def Copy(self):
        """Returns a copy of self."""
        return self.copy()
        
    def RawSource(self):
        """Returns raw source for this object."""
        s = ">>> inputmaps\n" + "\n".join(["%s = %s" % x for x in self.InputMaps]) + "\n"
        s += ">>> outputmaps\n" + "\n".join(["%s = %s" % x for x in self.OutputMaps]) + "\n"
        s += ">>> defs\n" + "\n".join(self.Defs) + "\n"
        s += ">>> externals\n" + "\n".join(self.Externals) + "\n"
        s += "\n".join([">>> %s\n%s" % (sn, s.RawSource()) 
                          for (sn,s) in self.Blocks.items()])
        return s
        
    def UpdateReserved(self):
        """Adds current input, output, defined, and external variables to the reserved table."""
        self.Reserved = [VarName(x) for (x,m) in self.InputMaps + self.OutputMaps]
        self.Reserved += [VarName(x) for x in self.Defs + self.Externals]
        self.Reserved = [x.lower() for x in self.Reserved]
        self.Reserved = [x for (i,x) in enumerate(self.Reserved)
                         if not x in self.Reserved[i+1:]]

    def CheckSource(self, Source, ErrOrigin = ""):
        """Checks that any tokens in s that are in Reserved are in all caps."""
        if len(Source) == 0: return
        for Token in AllTokens(Source):
            if Token.lower() in self.Reserved and Token.upper() != Token:
                s = "Token %s is reserved. Must put in all caps." % Token
                if len(ErrOrigin): s += "  Origin: %s." % ErrOrigin
                raise NameError(s)
                    
    def CheckDefs(self, Defs, ErrOrigin = ""):
        """Checks that no variables in Reserved are defined in Defs."""
        for v in Defs:
            v = VarName(v).lower()
            if v in self.Reserved:
                s = "Variable name %s is reserved and can't be defined." % v
                if len(ErrOrigin): s += "  Origin: %s." % ErrOrigin
                raise ValueError(s)        
                
    def CheckFortCode(self, FortCode, ErrOrigin = ""):
        """Checks another FortCode object for definitions and source capitalized."""
        self.CheckDefs(FortCode.Defs, ErrOrigin)
        for (BlockName, Block) in FortCode.Blocks.items():
            self.CheckSource(Block.RawSource(), ErrOrigin)
                        
    def AddVars(self, InputMaps = [], OutputMaps = [], Defs = [], Externals = [],
                SkipVarNames = [], Check = False):
        """Adds inputs/outputs/definitions/externals."""
        if Check: 
            self.CheckDefs(Defs)
        #add only unique vars
        self.InputMaps.extend([x for x in InputMaps if not x in self.InputMaps
                               and not VarName(x[0]) in SkipVarNames])
        self.OutputMaps.extend([x for x in OutputMaps if not x in self.OutputMaps
                               and not VarName(x[0]) in SkipVarNames])
        self.Defs.extend([x for x in Defs if not x in self.Defs 
                          and not VarName(x) in SkipVarNames])
        self.Externals.extend([x for x in Externals if not x in self.Externals 
                               and not VarName(x) in SkipVarNames])   

    def __AddLines(self, Block, lines, Dep, Indent, Check = False):
        """Parses lines to add."""
        if Block.lower() == "inputmaps":
            for line in lines:
                if "=" in line:
                    m = line.split("=")[:2]
                    m = (m[0].strip(), m[1].strip())
                    if not m in self.InputMaps:
                        self.InputMaps.append(m)
                elif len(line.strip()) > 0:
                    raise ValueError("Cannot parse InputMap %s" % line)
        elif Block.lower() == "outputmaps":
            for line in lines:
                if "=" in line:
                    m = line.split("=")[:2]
                    m = (m[0].strip(), m[1].strip())
                    if not m in self.OutputMaps:
                        self.OutputMaps.append(m)
                elif len(line.strip()) > 0:
                    raise ValueError("Cannot parse OutputMap %s" % line)
        elif Block.lower() == "defs":
            for line in lines:
                line = line.strip()
                if len(line) and not line in self.Defs:
                    self.Defs.append(line)
                    if Check:
                        self.CheckDefs([line])
        elif Block.lower() == "externals":
            for line in lines:
                line = line.strip()
                if len(line) and not line in self.Defs:
                    self.Externals.append(line)
        else:
            if not Block in self.Blocks:
                self.NewBlock(Block)
            sThis = "\n".join(lines)   
            self.Blocks[Block].Add(sThis, Dep, Indent)
            if Check:  
                self.CheckSource(sThis)

    def Add(self, Source = "", Dep = None, Indent = 0, Check = False, Block = None):
        """Adds additional code."""
        Source = CodeStrip(Source)
        if Block is None: 
            Block = self.CurBlock
        else:
            self.CurBlock = Block
        ThisLines = []
        for line in Source.split("\n"):
            if line.strip().startswith(">>>"):
                #last code
                self.__AddLines(Block, ThisLines, Dep, Indent, Check = Check)
                ThisLines = []
                #switch the source
                try:
                    Block = line.strip()[3:].split()[0]
                except (IndexError, ValueError):
                    raise ValueError("Could not process block identifier in:\n%s" % line)
            else:
                ThisLines.append(line)
        self.__AddLines(Block, ThisLines, Dep, Indent, Check = Check)
        
    def AddFortCode(self, other, Dep = None, Indent = 0, Check = False, 
                    Block = None, IfCondition = None):
        """Adds another FortCode to current.  Runs other.Process() first."""
        o = other.copy()
        o.Process()
        if IfCondition is None:
            Source = o.ProcessedSource
        else:
            Source = IfBracket(o.ProcessedSource, IfCondition)    
        if Block is None: Block = self.CurBlock
        self.Add(Source, Dep, Indent, Check, Block)
        self.AddVars(o.InputMaps, o.OutputMaps, o.Defs, o.Externals, Check=Check)
        
    def AddMissing(self, other, Check = False):
        """Adds blocks in other that are not in self."""
        for (BlockName, Block) in other.Blocks.items():
            if not BlockName in self.Blocks:
                self.Blocks[BlockName] = Block.copy()
        
    def AddFortCodeVars(self, other, SkipVarNames = [], Check = False):
        """Adds variables from other FortCode object."""
        o = other.Copy()
        o.ProcessVars()
        self.AddVars(o.InputMaps, o.OutputMaps, o.Defs, o.Externals, 
                     SkipVarNames = SkipVarNames, Check = Check)
   
    def Break(self, Block = None):
        """Adds a line feed to source."""
        if Block is None: Block = self.CurBlock
        self.Blocks[Block].Break()
           
    def AddCodeBlock(self, cb, Check = False, Block = None):
        """Adds other to current source code."""
        if Check:
            for (s, d, i) in cb.lines:
                self.CheckSource(s)
        if Block is None: Block = self.CurBlock
        self.Blocks[Block].AddCodeBlock(cb)
        
    def AddMapParent(self, Parent):
        """Adds a prefix to all mapped variables."""
        self.InputMaps = [(v, "%s.%s" % (Parent, m)) for (v,m) in self.InputMaps]
        self.OutputMaps = [(v, "%s.%s" % (Parent, m)) for (v,m) in self.OutputMaps] 
                 
    def ReplaceTokenSource(self, Token, Val):
        """Replaces Token with Val in sources."""
        for (BlockName, Block) in self.Blocks.items():
            Block.ReplaceTokenSource(Token, Val)
        if not self.ProcessedSource is None:
            self.ProcessedSource = ReplaceToken(self.ProcessedSource, Token, Val)

    def ReplaceTokenVars(self, Token, Val):
        """Replaces Token with Val in variable definitions."""
        self.InputMaps = [(ReplaceToken(v, Token, Val), ReplaceToken(m, Token, Val))
                          for (v,m) in self.InputMaps]
        self.OutputMaps = [(ReplaceToken(v, Token, Val), ReplaceToken(m, Token, Val))
                           for (v,m) in self.OutputMaps]
        self.Defs = [ReplaceToken(v, Token, Val) for v in self.Defs]
        self.Externals = [ReplaceToken(v, Token, Val) for v in self.Externals]                 

    def __GetHolders(self, BlockName, Prev = []):
        """Returns all holders that BlockName depends on."""
        HolderDep = set()
        for bn in self.Blocks:
            if bn == BlockName: continue
            if self.HasHolder[BlockName][bn]:
                if bn in Prev:
                    raise ValueError("Found circular dependencies for %s" % bn)
                else:
                    HolderDep.add(bn)
                    HolderDep.update(self.__GetHolders(bn, Prev + [BlockName]))
        return HolderDep
        
    def ProcessSource(self):
        """Processes and returns source code."""
        #maks a list of which holders each source has
        self.HasHolder = {}
        for bn1 in self.Blocks:
            self.HasHolder[bn1] = {}
            for bn2 in self.Blocks:
                if bn1==bn2: continue
                self.HasHolder[bn1][bn2] = self.Blocks[bn1].HasHolder(bn2)
        #now get dependencies
        self.HolderDep = {}
        for bn1 in self.Blocks:
            self.HolderDep[bn1] = self.__GetHolders(bn1)
        #sort dependencies by least dependent first
        Dep = [(len(v), k, v) for (k,v) in self.HolderDep.items()]
        Dep.sort()
        #now process each source in series and build up code
        HolderDict = {}
        for (lenv, bn, v) in Dep:
            cb = self.Blocks[bn].Copy()
            HolderDict[bn] = cb.Process(HolderDict = HolderDict)
        #return the main code
        self.ProcessedSource = HolderDict["Main"]
        #check if condition
        if not self.IfCondition is None and len(self.ProcessedSource.strip()):
            self.ProcessedSource = IfBracket(self.ProcessedSource, self.IfCondition)
        #add comment
        if len(self.ProcessedSource.strip()) and len(self.Comment):
            self.ProcessedSource = "!%s\n" % self.Comment + self.ProcessedSource
        return self.ProcessedSource
        
    def ProcessMain(self, *args):
        """Returns processed source with Main replaced with source from args."""
        o = self.copy()
        o.Blocks["Main"].Clear()
        for itm in args:
            o.SetCurrent("Main")
            o.Add(itm)
        return o.Process()  
  
    def NewMain(self, *args):
        """Returns a copy of self with Main replaced with source from args."""
        o = self.copy()
        o.Blocks["Main"].Clear()
        for itm in args:
            o.SetCurrent("Main")
            o.Add(itm)
        return o  
       
    def ProcessVars(self):
        """Adds an index to variables that have two underscores after variable name."""
        Vars =   [v for (v,m) in self.InputMaps] \
               + [v for (v,m) in self.OutputMaps] \
               + [v for v in self.Defs]
        for Var in Vars:
            CBlock, Type, Name, Shape, Val = VarSplit(Var)
            if Name.endswith("__"):
                if self.VarIndex is None:
                    raise ValueError("Index cannot be added to variable %s" % Name)
                else:
                    NewName = Name[:-1] + "%d" % self.VarIndex
                    self.ReplaceTokenSource(Name, NewName)
                    self.ReplaceTokenVars(Name, NewName)
         
    def Process(self):
        """Removes duplicate definitions, etc, from all variable lists.
Also checks for any multiply-defined maps.  Also removes variables
that are not used in source."""
        #functions to check all of the variable definitions
        def Prune(l):
            return [x for (i,x) in enumerate(l) if not x in l[i+1:]]
        def CheckMapDuplicates(l):
            VarDefs = [v for (v,m) in l]
            for Var in VarDefs:
                if VarDefs.count(Var) > 1:
                    s = "\n".join(["%s = %s" % (v,m) for (v,m) in l if v == Var])
                    raise ValueError("Variable %s is mapped in more than one way:\n%s" % (Var,s))
        #first process the source
        self.ProcessSource()
        #process variables
        self.ProcessVars()
        #move any inout variables to the end of Inputs
        OutputNames = [VarName(Output) for (Output, Map) in self.OutputMaps]
        InOutMaps = [(Input, Map) for (Input, Map) in self.InputMaps
                     if VarName(Input) in OutputNames]
        self.InputMaps = [InputMap for InputMap in self.InputMaps
                          if not InputMap in InOutMaps] + InOutMaps        
        #prune and check for duplicates
        self.InputMaps = Prune(self.InputMaps)
        CheckMapDuplicates(self.InputMaps)
        self.OutputMaps = Prune(self.OutputMaps)
        CheckMapDuplicates(self.OutputMaps)
        self.Defs = Prune(self.Defs)
        self.Externals = Prune(self.Externals)
        #make a dependency tree
        Depend = {}
        Keep = []
        AllVars = [v for (v,m) in self.InputMaps + self.OutputMaps] + self.Defs
        AllVars = [v.lower() for v in AllVars]
        for v in AllVars:
            vn = VarName(v)
            if not vn in Depend:
                Depend[vn] = []
            for vs in VarShape(v):
                vsn = VarName(vs)
                if not vsn in Depend:
                    Depend[vsn] = []
                if not vsn in Depend[vn]:
                    Depend[vn].append(vsn)
        for vn in Depend.keys():
            if HasToken(self.ProcessedSource, vn):
                Keep.append(vn)
        #find variables to keep
        LastLen = -1
        while LastLen != len(Keep):
            LastLen = len(Keep)
            for vn in Keep:
                for vsn in Depend[vn]:
                    if not vsn in Keep:
                        Keep.append(vsn)
        #get rid of variables that are not in source
        def Test(Var):
            return VarName(Var).lower() in Keep
        self.InputMaps = [(v,m) for (v,m) in self.InputMaps
                          if Test(v)]
        self.OutputMaps = [(v,m) for (v,m) in self.OutputMaps
                           if Test(v)]
        self.Defs = [v for v in self.Defs if Test(v)]
        #return value
        return self.ProcessedSource
        

TEST = False
if TEST:
    s = """
a = 1
b = 2
c = 4      || [sub1]
>>> sub1
h = 5
>>> Main
d = 6      || k
f = 7      || l
[sub1]
[sub2]
>>> sub2
k = 3
>>> sub3
l = 5"""
    fc = FortCode(Source = s)
    for (k,v) in fc.Blocks.items():
        print k, "\n", v.Lines
    print "=="
    print fc.Process()
    exit()

