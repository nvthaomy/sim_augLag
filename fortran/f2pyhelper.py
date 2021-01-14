#/usr/bin/env python

### Helper routines for compiling with f2py

import os, sys


#check if there is a configuration file
Args = ""
if "F2PYARGS" in os.environ:
    Args = os.environ["F2PYARGS"]
elif os.path.isfile("f2py.arg"):
    Args = file("f2py.arg").readline().split()
else:
    Paths = [os.path.dirname(__file__)] + sys.path 
    for Dir in Paths:
        File = os.path.join(Dir, "sim/f2py.arg")
        if os.path.isfile(File):
            Args = file(File).readline().split()
            break
        File = os.path.join(Dir, "f2py.arg")
        if os.path.isfile(File):
            Args = file(File).readline().split()
            break

#Should NOT use InPlace because slices will point to unallocated space after passing to fortran routines 
UseInplace = False
if "--inplace" in Args:
    UseInplace = True
    Args = [x for x in Args if not x == "--inplace"]
if "--noinplace" in Args:
    UseInplace = False
    Args = [x for x in Args if not x == "--noinplace"]
 
ReportOnCopy = -1
for (i, Arg) in enumerate(Args):
    if Arg.startswith("--reportoncopy="):
        ReportOnCopy = int(Arg[15:])




def compile(source,
            modulename = 'untitled',
            extra_args = None,
            source_fn = None
            ):
    ''' Build extension module from processing source with f2py.
    Read the source of this function for more information.
    '''
    import subprocess
    import tempfile
    import os, sys
    
    if extra_args is None:
        extra_args = Args
    if type(extra_args) is list:
        extra_args = " ".join(extra_args)
        
    if source_fn is None:
        fname = os.path.join(tempfile.mktemp()+'.f')
    else:
        fname = source_fn

    f = open(fname,'w')
    f.write(source)
    f.close()
    
    Pythonexe = sys.executable

    args = ' -c -m %s %s %s'%(modulename,fname,extra_args)
    if ReportOnCopy > 0:
        args = args + " -DF2PY_REPORT_ON_ARRAY_COPY=%d" % ReportOnCopy
    c = '%s -c "import numpy.f2py as f2py2e;f2py2e.main()" %s' %(Pythonexe,args)

    p = subprocess.Popen(c, shell = True,
                         stdout = subprocess.PIPE,
                         stderr = subprocess.STDOUT)
    o = p.stdout.read()
    p.stdout.close()
    p.wait()
    s = p.returncode

    if source_fn is None:
        try: os.remove(fname)
        except OSError: pass
    return s, o                     
    