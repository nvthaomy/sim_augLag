#/usr/bin/env python

### Imports the SIM package.

import numpy as np


#these need to be imported in this order
#to address dependencies 
import units
import utility
import histogram
import mathutil
import spline
import fortran

import geom
import chem
import potential

import measure
import atomselect
import atommap
import integrate
import traj

import system
import srel
import cluster

import alg
import export


#import a settings file if given; search all paths including current
import sys, os
Paths = [os.getcwd()] + sys.path + [os.path.dirname(__file__)]
for p in Paths:
    fn = os.path.join(p, "simsettings.py")
    if os.path.isfile(fn):
        execfile(fn)
        break
del sys
del os
