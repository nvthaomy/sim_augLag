#/usr/bin/env python

### Fortran module generation


import f2pyhelper


MaxLineLen = 130
            
        
class Module(list):
    
    def __init__(self, Name, Path = None, ForceCompile = False, KeepSource = False, 
                 UseFortModule = False, FortModuleName = "m"):
        """Defines a new module of subroutines in Fortran.
"""
        import os
        object.__setattr__(self, "_Module__Vars", {}) 
        self.__Vars0 = []
        list.__init__(self)
        self.Code = ""
        self.Path = Path
        if self.Path is None:
            self.Path = "./pylib"
            if not os.path.isdir(self.Path):
                os.mkdir(self.Path)
        self.KeepSource = KeepSource
        self.ForceCompile = ForceCompile
        #names
        self.Name = Name
        #pre and post loading functions
        self.PreLoadList = []
        self.PostLoadList = []
        #compiled module
        self.Module = None
        self.BaseModule = None
        #fortran module structure
        self.UseFortModule = UseFortModule
        self.FortModuleName = FortModuleName
        #source
        self.Source = None

        
    def __getattr__(self, name):
        """Allows direct accessing of variables."""
        Var, IsConstant = self.__Vars.get(name, (None,None))
        if Var:
            if not Var.shape: Var = Var.item()
            return Var
        else:
            raise AttributeError(name)
            
    def __setattr__(self, name, value):
        """Allows direct accessing of variables."""
        Var, IsConstant = self.__Vars.get(name, (None,None))
        if Var:
            if IsConstant:
                raise AttributeError("%s is a constant and cannot be set." % name)
            elif Var.shape:
                Var[:] = value
            else:
                Var.itemset(value)
        else:
            object.__setattr__(self, name, value)

        
    def AddVars(self, VarStr, Prefix = "", IsConstant = False):
        """Adds module-level variables."""
        from var import VarSplit
        VarList = []
        VarStr = VarStr.strip()
        for VarSpec in VarStr.split("\n"):
            if not VarSpec: continue
            CBlock, Type, Name, Shape, Val = VarSplit(VarSpec)
            if any([x[0]==Name for x in self.__Vars0]):
                raise AttributeError("%s already defined for Lib object." % Name)
            Name = Prefix + Name 
            self.__Vars0.append((Name, VarSpec, IsConstant))
            VarList.append(Name)
        return VarList
            
    def AddVarFromObj(self, VarName, Obj, Prefix = "", IsConstant = False, ShapeStr = None):
        """Adds a module-level constant from a Python object."""
        from var import GetVarSpec, VarSplit
        spec = GetVarSpec(VarName, Obj, IsConstant = IsConstant, ShapeStr = ShapeStr)
        CBlock, Type, Name, Shape, Val = VarSplit(spec)
        Name = Prefix + Name 
        if any([x[0]==Name for x in self.__Vars0]):
            raise AttributeError("%s already defined for Lib object." % Name)        
        self.__Vars0.append((Name, spec, IsConstant))       
        return Name
    
    
    def VarPreLoad(self, Obj, Vars = [], Consts = [], Prefix = ""):
        """Provess variables before compilation."""
        if not hasattr(Obj, "LibVars0"):
            Obj.LibVars0 = {}         
        #aggregate all of the variables
        AllVars = [(x, True) for x in Consts] + [(x, False) for x in Vars]
        #constants from various objects
        for (itm, IsConstant) in AllVars:
            #check for a tuple
            if isinstance(itm, (tuple,list)):
                if not len(itm) == 2:
                    raise ValueError("Expected 2 arguments from tuple %s" % str(itm))
                ObjName, Name = itm
            else:
                ObjName = itm
                Name = itm
            ret = self.AddVarFromObj(Name, getattr(Obj, ObjName), 
                                     Prefix = Prefix, IsConstant=IsConstant)
            Obj.LibVars0[ObjName] = ret

    def VarPostLoad(self, Obj):
        """Run after compilation."""
        HasLibVars = hasattr(Obj, "LibVars")
        for (ObjName, ModuleName) in Obj.LibVars0.items():
            if hasattr(Obj, ObjName):
                #skip constants
                ModuleVar, IsConstant = self.__Vars.get(ModuleName)
                if IsConstant: continue
                #copy the value of the variable to the library version
                ObjVar = getattr(Obj, ObjName)
                #now delete the old variable
                delattr(Obj, ObjName)
                #initialize values
                if ModuleVar.shape:
                    ModuleVar[:] = ObjVar
                else:
                    ModuleVar.itemset(ObjVar)
                #initialize the new variable and make a link
                if HasLibVars:   
                    Obj.LibVars[ObjName] = ModuleVar
                setattr(Obj, ObjName, ModuleVar)
        #update so we automatically get attributes by grabbing the lib ones
        del Obj.LibVars0
        
        
       
    def GetFortran(self, Force = False):
        """Makes the Fortran code for this module."""
        from var import VarSplit, IsAllocatable
        from code import I
        #check if we've already made the source
        if not self.Source is None and not Force:
            return self.Source
        #reorder variables to have constants first
        self.__Vars0 = [x for x in self.__Vars0 if x[-1]] + [x for x in self.__Vars0 if not x[-1]]
        DoneSubs = {}
        s = "!======MODULE CODE FOR %s======\n\n" % self.Name
        if self.UseFortModule:
            s += "module %s\n\n" % self.FortModuleName
            #do the module-level constants
            for (Name, Var, IsConstant) in self.__Vars0:
                CBlock, Type, OldName, Shape, Val = VarSplit(Var, FortranType = True)
                if len(CBlock):
                    raise ValueError("Error in specification for module level variable in a common block.")
                t = "%s" % Type
                if IsAllocatable(Shape):
                    if IsConstant:
                        raise ValueError("Constant specification cannot be allocatable.")
                    else:
                        t += ", allocatable"
                if Shape:
                    t += ", dimension(%s)" % Shape
                if IsConstant:
                    t += ", parameter :: %s" % Name
                else:
                    t += " :: %s" % Name
                if Val:
                    t += " = %s\n" % Val
                s += I(t, 4, MaxLen = MaxLineLen)                    
            s += "\n\ncontains\n\n"
        else:
            if len(self.__Vars0):
                raise ValueError("Found Fortran module-level variables but UseFortModule = False.")
        for Sub in self:
            if Sub.Name in DoneSubs:
                raise NameError("Found multiple definitions for %s Fortran routine." % Sub.Name)
            for DoneSub in DoneSubs.values():
                Sub.UpdateExternal(DoneSub)
            Sub.Process()
            s += I(Sub.ProcessedSource, MaxLen = MaxLineLen)
            s += "\n\n"
            #add to done subs list
            DoneSubs[Sub.Name] = Sub
        s += self.Code
        if self.UseFortModule:
            s += "\n\nend module\n"
            #get rid of all of the "external" references because they will be accessible
            #automatically within the fortran module structure
            news = ""
            for line in s.split("\n"):
                if not line.strip().lower().startswith("external "):
                    news += line + "\n"
            s = news
        return s


    def Hash(self):
        """Returns a hash string."""
        import hashlib
        return hashlib.sha1(self.GetFortran()).hexdigest()


    def Unload(self):
        """Unloads any modules or compiled functions from memory."""
        import sys
        if not self.Module is None:
            self.__Vars.clear()
            if self.Name in sys.modules:
                del sys.modules[self.Name]
            for Sub in self:
                Sub.Func = None
            if not self.BaseModule is None: 
                del self.Module
                self.Module = self.BaseModule
                del self.BaseModule
            nrefs = sys.getrefcount(self.Module)
            if nrefs > 2:
                print "WARNING: cannot unload module %s -> too many references (%d)." % (self.Name, nrefs)
            self.Module = None


    def __GetModuleFile(self):
        """Returns the name of the compiled module file."""
        import imp, os
        self.__ModuleFile = None
        for (Suffix, Mode, Type) in imp.get_suffixes():
            if Type == imp.C_EXTENSION:
                fn = os.path.join(self.Path, "%s%s" % (self.Name, Suffix))
                if os.path.isfile(fn):
                    self.__ModuleFile = os.path.normpath(fn)
                    break


    def __PreLoad(self):
        """Runs any pre-load functions."""
        for Sub in self:
            Sub.PreLoad()  
        for PreLoad in self.PreLoadList:
            PreLoad(self)


    def __PostLoad(self):
        """Links subroutine functions to the module and runs any post-load functions."""
        for Sub in self:
            Sub.Func = getattr(self.Module, Sub.Name)
            setattr(Sub, Sub.Name, Sub.Func)
            SubStr = Sub.SubStr
            if len(SubStr):
                SubStr = SubStr + "."
            #make a function call string
            if len(Sub.OutputStr):
                Sub.CallStr = "%s = %s%s(%s)" % (Sub.OutputStr, SubStr,
                                                 Sub.Name, Sub.InputStr)
            else:
                Sub.CallStr = "%s%s(%s)" % (SubStr, Sub.Name, Sub.InputStr)                
            if len(Sub.PreCallCode):
                Sub.CallStr = Sub.PreCallCode + "\n" + Sub.CallStr
            Sub.CallObj = compile(Sub.CallStr, '<string>', 'exec') 
        #run all post-load functions
        for Sub in self:
            Sub.PostLoad()   
        for PostLoad in self.PostLoadList:
            PostLoad(self)


    def Load(self, Verbose = 1):
        """Loads the current module."""
        self.__GetModuleFile()
        self.__PreLoad()
        if self.__ModuleFile is None or self.ForceCompile:
            self.__Compile()
        else:
            import imp
            if Verbose:
                print "> Loading compiled library %s" % self.__ModuleFile
            #try to load existing library
            try:
                self.Module = imp.load_dynamic(self.Name, self.__ModuleFile)
            except ImportError:
                print "> Found existing library but could not import."
                raise ImportError
            else:
                UpToDate = False
                try:
                    sHash = self.Module.modulehash()
                except (AttributeError, NameError):
                    pass
                else:
                    UpToDate = (sHash == self.Hash())
                if not UpToDate:
                    s = """Module %s is out of date!
Delete compiled file %s and re-run script to automatically recompile.""" % (self.__ModuleFile, self.Name)
                    raise ImportError(s)
        #transforme Vars data for faster use of __getattr__ and __setattr__
        self.__Vars = {}
        for (Name, Var, IsConstant) in self.__Vars0:
            Var = getattr(self.Module, Name.lower())
            self.__Vars[Name] = (Var, IsConstant)
        self.__PostLoad()
              
            
    def __Compile(self, Verbose = 1):
        """Compiles a Fortran module."""
        import os, imp
        
        #add a hash function
        sHash = self.Hash()
        self.Code = self.Code + """

function modulehash()
    character(len=%d) :: modulehash
    modulehash = '%s'
end function
""" % (len(sHash), sHash)

        #get source
        self.Source = self.GetFortran()
                   
        #make a filename and compile
        SrcFile = os.path.join(self.Path, self.Name + ".f90")
        SrcFile = os.path.abspath(SrcFile)

        if Verbose:
            print "> Compiling Fortran for %s " % self.Name

        cwd = os.getcwd()
        os.chdir(self.Path)
        Status, Output = f2pyhelper.compile(self.Source,
                                            modulename = self.Name,
                                            source_fn = SrcFile)
        os.chdir(cwd)
        
        #check status
        if Status > 0:
            print Output
            print "!"*20
            s = "Could not compile %s" % SrcFile
            s += "\nModule file may already be in use."
            raise StandardError(s)
        if Verbose == 2:
            print Output
        if Verbose:
            print "> Successfully compiled %s " % self.Name
        if not self.KeepSource:
            os.remove(SrcFile)

        #load the module object
        self.__GetModuleFile()
        if self.__ModuleFile is None:
            raise StandardError("Could not find compiled module for %s" % self.Name)
        else:
            try:
                self.Module = imp.load_dynamic(self.Name, self.__ModuleFile)
            except ImportError, Val:
                print Val
                s = "Could not load compiled module %s for %s" % (self.__ModuleFile, self.Name)
                raise StandardError(s)
            if self.UseFortModule:
                self.BaseModule = self.Module
                self.Module = getattr(self.Module, self.FortModuleName)

        
          
    

            


        