#/usr/bin/env python

### Initialization file for measuring components

import glob, os

import base
from base import MeasuresClass, ListMeasures
import standard
from standard import *


Query = os.path.join(__path__[0], "*.py")

for fn in glob.glob(Query):   
    
    Package, Ext = os.path.splitext(os.path.basename(fn))
    if not Package in ["__init__", "base", "basefortran", "standard"]:
        s = "from %s import *" % Package
        exec(s)