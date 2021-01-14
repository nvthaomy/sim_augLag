# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 14:56:37 2017

@author: Shell
"""


import sim
import numpy as np

Lib = sim.fortran.Module("testmodule", UseFortModule = True, ForceCompile = True, 
                         KeepSource = True, Path = "./")



Lib.Code = """
subroutine test0()
    print *, "TEST0"
    print *, a
    print *, b
    print *, c
    print *, d
    print *, lbound(b), ubound(b)
    print *, lbound(c), ubound(c)
    print *, size(pos,1), size(pos,2)
    print *, lbound(pos,1), ubound(pos,1)
    print *, lbound(pos,2), ubound(pos,2)
    print *, "TEST0 done."
end subroutine

subroutine test1()
    print *, "TEST1"
    pos = 0.d0
end subroutine

subroutine test2()
    print *, "TEST2"
    print *, sum(pos, dim=1) / size(pos, dim=1)
end subroutine
"""
Lib.AddVars("""
float pos(:,:)
int n = 5""")
Lib.AddVars("""float a = 5.
float b(0:2) = (/ 3, 4, 9 /)
bool c(1:4) = (/ 0, 1, 1, 0 /)
""", IsConstant = True)
d = np.arange(6).reshape((2,3))
Lib.AddVarFromObj("d", d, IsConstant = True)
Lib.Load()


print type(Lib.Module.n)
print Lib.Module.n.dtype
n2 = Lib.Module.n
print n2
Lib.Module.n = 10
print Lib.Module.n, n2

print Lib.n
print type(Lib.n)
    

Lib.Module.pos = [[1,2,3],[4,5,6]]  #np.arange(6).reshape((2,3))
Lib.Module.test0()
Lib.Module.pos = np.arange(10).reshape((5,2))
Lib.Module.test0()
print Lib.Module.pos.shape

print Lib.BaseModule.__doc__
print Lib.Module.__doc__

import sys
print sys.getrefcount(Lib.Module)
Lib.Module.test0()
print "Python values"
print Lib.Module.a
print Lib.Module.b
print Lib.Module.c
print Lib.Module.d

Lib.Module.test1()

Lib.Module.test2()

Lib.Module.pos[:] = np.arange(6).reshape((2,3))

print Lib.Module.pos

Lib.Module.test2()

Pos = Lib.Module.pos
Pos += 1.
print Pos
Pos = Pos + 1.
print Lib.Module.pos

print sys.getrefcount(Lib.Module)
del Pos
print sys.getrefcount(Lib.Module)



