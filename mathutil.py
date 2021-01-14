#!/usr/bin/env python


import numpy as np


def LeastSquaresFit(Xs, Ys, Weights = None):
    """Fits a linear model to a curve of interest."""
    n, m = Xs.shape
    if not len(Ys) == n:
        raise TypeError("Size of Ys does not match XCoefs")
    if Weights is None:
        Weights = np.ones(n, float)
    elif not len(Weights) == n:
        raise TypeError("Size of Weights does not match XCoefs")
    W = np.diag(np.sqrt(Weights))
    XW = np.dot(W, Xs)
    YW = np.dot(Ys, W)
    Coefs = np.linalg.lstsq(XW, YW)[0]
    return Coefs
    
   


