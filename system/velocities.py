#/usr/bin/env python


### Routines for initial velocities in SIM suite.
### coded by MSS

import numpy as np


def Canonical(Sys, Temp):
    """Gives Maxwell-Boltzmann velocities to system, with zero COM."""
    Sigma = np.sqrt(Temp * Sys.Units.kB) * Sys.sqrtMass 
    p = np.random.normal(size = Sys.Vel.shape) * Sigma[:,np.newaxis]
    p = p - p.mean(axis=0)
    Sys.Vel = p * Sys.iMass[:,np.newaxis]
    
    
    
        
