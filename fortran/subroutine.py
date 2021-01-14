#/usr/bin/env python

### Fortran subroutine generation


from var import VarShape, VarNorm, VarSplit, VarName
from code import Unique, I, MaxLen, HasToken

import f2pyhelper



class Subroutine(list):
    
    def __init__(self, Name, FortCode = None, SubStr = "", PreCallCode = ""):
        """Defines a new subroutine.
Arguments:
    Name: name of subroutine
    FortCode: a FortCode object 
    SubStr: string name of this sub variable"""
        list.__init__(self)
        self.Name = Name.lower()
        #reference to the name of this variable
        self.SubStr = SubStr   
        #pre commands to calling
        self.PreCallCode = PreCallCode
        #processed variables for fortran
        self.AddFortCode(FortCode)
        #strings for calling in Python
        self.InputStr = None
        self.OutputStr = None
        #string for calling in Fortran
        self.FortCallStr = None
        #processed source
        self.ProcessedSource = None
        #compiled function
        self.Func = None        
        
    def AddFortCode(self, FortCode):
        """Adds FortCode object."""
        if FortCode is None:
            self.FortCode = None
        else:
            self.FortCode = FortCode.copy()
    
    
    def UpdateExternal(self, RefSub):
        """If this sub references an external with [external:functionname], 
then automatically update it with the relevant information."""  
        Block = "external:%s" % RefSub.Name
        CallStr = "call %s" % RefSub.FortCallStr
        fc = self.FortCode
        if self.FortCode.ReferencesBlock(Block):
            #delete the old block if it exists
            if fc.HasBlock(Block):
                fc.DelBlock(Block)
            #add only vars that are not in Defs
            SkipVarNames = [VarName(x) for x in fc.Defs]
            fc.AddFortCodeVars(RefSub.FortCode, SkipVarNames = SkipVarNames)
            fc.Add(CallStr, Block = Block)
            fc.AddVars(Externals = [RefSub.Name])
            
    def Process(self):
        """Makes the Fortran code for this subroutine."""
        #process the fortran code
        self.Source = self.FortCode.Process()
        self.InputMaps = self.FortCode.InputMaps
        self.OutputMaps = self.FortCode.OutputMaps
        self.Defs = self.FortCode.Defs
        self.Externals = self.FortCode.Externals
                      
        #process the maps
        self.Inputs = [Input for (Input, Map) in self.InputMaps]
        self.Outputs = [Output for (Output, Map) in self.OutputMaps]

        #add dimension variables to inputs
        Shapes = []
        for x in self.Inputs:
            Shapes.extend(VarShape(x))
            
        #remove any shape variables that are in inputs
        ShapeNames = Unique([VarName(x) for x in Shapes])
        self.InputMaps = [(Input, Map) for (Input,Map) in self.InputMaps 
                          if not VarName(Input) in ShapeNames]
        self.Inputs = [Input for (Input, Map) in self.InputMaps]
        
        #check that no shape variables in output
        for x in self.Outputs:
            if VarName(x) in Shapes:
                raise ValueError("%s is a shape variable -- cannot be in output variable list." % x)

        #normalize
        Inputs = [VarNorm(x) for x in self.Inputs]
        Outputs = [VarNorm(x) for x in self.Outputs]
        Defs = [VarNorm(x) for x in self.Defs]

        #arg line      
        Args = Unique(Inputs + Outputs + Shapes)       
        #defs
        Defs = Unique(Shapes + Inputs + Outputs + Defs)
        
        #sort arguments
        Inputs.sort(key = lambda x: Args.index(x)) 
        Outputs.sort(key = lambda x: Args.index(x)) 

        #start the code
        ArgLine = ", &\n     & ".join([VarName(x) for x in Args])
        self.FortCallStr = "%s(%s)" % (self.Name, ArgLine)
        s = I("subroutine %s(%s)\n" % (self.Name, ArgLine))
        s += I("implicit none", 4)
        sCommon = ""

        for Def in Defs:
            CBlock, Type, Name, Shape, Val = VarSplit(Def, FortranType = True)
            sShape = ""
            if len(Shape):
                Shape = Shape.replace(" ","")
                if Shape.startswith("1:"):
                    Shape = Shape[2:]
                sShape = ", dimension(%s)" % Shape
            if len(CBlock):
                sCommon += "COMMON /%s/ %s\n" % (CBlock, Name)
            sIntent, sVal = "", ""
            InOut = False
            Out = False
            if Def in Shapes:
                sIntent = ", intent(in)"
            elif Def in Inputs:
                if Def in Outputs:
                    sIntent = ", intent(inout)"
                    InOut = True
                else:
                    sIntent = ", intent(in)"
            elif Def in Outputs:
                sIntent = ", intent(out)"
                Out = True
            elif len(Val):
                sIntent = ", parameter"
                sVal = " = %s" % Val
            s += I("%s%s%s :: %s%s" % (Type, sShape, sIntent, Name, sVal), 4)
            if InOut:
                if len(Shape) > 1 and f2pyhelper.UseInplace:
                    s += I("!f2py intent(in, out, inplace) :: %s" % Name, 0)
                else:
                    s += I("!f2py intent(in, out) :: %s" % Name, 0)
            if Out:
                s += I("!f2py intent(out) :: %s" % Name, 0)   
                
        #add common block
        if len(sCommon):
            s += I(sCommon, 4)                 
        #add externals
        for Ext in self.Externals:
            if HasToken(self.Source, Ext.strip()):
                s += I("external :: %s" % Ext, 4)

        #add source code
        s += "\n"
        s += I(self.Source, 4)
        s += "\n"
        s += I("end subroutine")

        #set input and output strings -- processing here is to put in the same order as the Arg line
        InputMaps = dict([(VarNorm(Input), Map) for (Input,Map) in self.InputMaps])
        OutputMaps = dict([(VarNorm(Output), Map) for (Output,Map) in self.OutputMaps])
        self.InputStr = ", ".join([InputMaps[Input] for Input in Inputs])
        self.OutputStr = ", ".join([OutputMaps[Output] for Output in Outputs])
        
        self.ProcessedSource = I(s, MaxLen = MaxLen)
        
        return self.ProcessedSource
   


    def PreLoad(self):
        """Function that is called before compiling."""
        pass

    def PostLoad(self):
        """Function that is called after compiling."""
        pass
    