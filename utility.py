#!/usr/bin/env python

import sys
import time
import tempfile
import os
import glob


ProgressBarOn = sys.stdout.isatty()


def NewTempFile(Suffix = '', Dir = None):
    """Returns the name of a new temporary file."""
    Fd, Fn = tempfile.mkstemp(suffix=Suffix, dir=Dir)
    os.close(Fd)
    return Fn
    
def NewTempFilePrefix(Dir = None):
    """Returns the name of a new temporary file with name prefix"""
    Fd, Fn = tempfile.mkstemp(suffix="", dir=Dir)
    os.close(Fd)
    os.remove(Fn)
    return Fn    


class ProgressBar(object):
    def __init__(self, Text, Steps = 1, BarLen = 20, UpdateFreq = 1.):
        """Initializes a generic progress bar."""
        self.Text = Text
        self.Steps = Steps
        self.BarLen = BarLen
        self.UpdateFreq = UpdateFreq
        self.__LastTime = 0.
        self.__LastTime = time.time()
        self.__LastLen = 0
        self.Update(0)

    def Update(self, Step):
        """Updates the progress bar."""
        if time.time() - self.__LastTime > self.UpdateFreq:
            if not ProgressBarOn:
                return
            self.__LastTime = time.time()
            if self.BarLen == 0:
                s = "%s [%d]" % (self.Text, Step)
            else:
                Frac = float(Step) / (self.Steps + 1.e-300)
                n = int(self.BarLen * Frac + 0.5)
                n = max(min(n, self.BarLen), 0)
                s = "%s [" % self.Text
                s += "="*n + (self.BarLen-n)*" "
                s += "] %.1f%%" % (100.*Frac)
            self.__LastLen = len(s)
            s += "\r"
            sys.stdout.write(s)
            sys.stdout.flush()

    def Clear(self):
        """Clears text on this line."""
        if not ProgressBarOn:
            return
        sys.stdout.write(" "*self.__LastLen + "\r")
        sys.stdout.flush()
       
       
       
       
       
def ExpandArgs(Args, NTOnly = False):
    """Expands arguments that contains wildcards.
On Linux systems, does nothing since command interpreter
already does this."""
    if os.name == "nt" or not NTOnly:
        NewArgs = []
        for a in Args:
            if ("*" in a or "?" in a) and not a.startswith("--"):
                NewArgs.extend(glob.glob(a))
            else:
                NewArgs.append(a)
        return NewArgs
    else:
        return Args

def ArgAsType(Val, Template):
    """Returns Val as the type of Template"""
    #for none we will return True
    if Template is None:
        return True
    else:
        t = type(Template)
        return t(Val)
      
def ParseArgs(Defaults = {}, Warn = False, AcceptArg = None):
    """Parses a list of arguments into a dictionary.
The keys 0,1,..n give the consecutive non-option arguments.
The key "FLAGS" gives a list of flags used.
The key "FLAGPOS" tells the number of non-flag arguments before each flag.
The key "ARGS" gives a list of non-flag arguments.
The key "NARG" gives the number of arguments.
Other keys give the values of flag arguments.
Flags can be '-x', '--word', or '--word=value'."""
    ind = 0
    Ret = Defaults.copy()
    Ret["FLAGS"] = []
    Ret["FLAGPOS"] = []
    Argv = ExpandArgs(sys.argv)
    if AcceptArg is None:
        AcceptArg = []
    for k in Defaults.keys():
        if not k in AcceptArg: AcceptArg.append(k)
    for i in range(0, len(Argv)):
        a = Argv[i]
        if a.startswith("--"):
            #word flags
            if "=" in a:
                Flag = a[2:].split("=")[0]
                Val = a[2:].split("=")[1]
                HasVal = True
            else:
                #get the following arg
                Flag = a[2:]
                if i + 1 < len(Argv):
                    Val = Argv[i+1]
                else:
                    Val = ""
                HasVal = False
            #try to get the same type as defaults
            if Flag in Defaults:
                Ret[Flag] = ArgAsType(Val, Defaults[Flag])
            else:
                Ret[Flag] = Val
            if Warn and HasVal and not Flag in AcceptArg:
                print "Did not recognize option %s" % Flag
            Ret["FLAGS"].append(Flag)
            Ret["FLAGPOS"].append(ind)
        elif a.startswith("-") and len(a) > 1 and not a[1] in "0123456789":
            #one letter flags, use following arg as value
            if i + 1 < len(Argv):
                Val = Argv[i+1]
            else:
                Val = ""
            for Flag in a[1:]:
                #try to get the same type as defaults
                if Flag in Defaults:
                    Ret[Flag] = ArgAsType(Val, Defaults[Flag])
                else:
                    Ret[Flag] = Val
                Ret["FLAGS"].append(Flag)
                Ret["FLAGPOS"].append(ind)
        else:
          #normal argument
          Ret[ind] = a
          ind += 1
    Ret["NARG"] = ind
    Ret["ARGS"] = [Ret[i] for i in range(ind)]
    return Ret


def ParseNumList(Arg, Offset = 0):
    """Allows number specifications like '1-2,5,60-70'.
Returns a list, or None if Arg does not match anything."""
    if Arg is None or Arg == "": return None
    Ret = []
    for s in """'"[]""":
        Arg = Arg.replace(s,"")
    for l in Arg.split(","):
        if "-" in l:
            a, b = [int(x) for x in l.split("-")]
            Ret.extend(range(a + Offset, b+1+Offset))
        else:
            a = int(l)
            Ret.append(a+Offset)
    Ret.sort()
    Ret = [x for (i,x) in enumerate(Ret) if not x in Ret[i+1:]]
    return Ret