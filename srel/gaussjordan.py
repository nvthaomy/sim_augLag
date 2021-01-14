#/usr/bin/env python


### Gauss-Jordan elimination.
### coded by MSS

import numpy as np


def ReduceMatrix(M, Tol = 1.e-200):
    """Removes all-zero rows from a matrix."""
    M2 = []
    for Row in M:
        if not np.all(np.abs(Row[:-1]) < Tol):
            M2.append(Row)
    return np.array(M2, float)


def GaussJordan(M, Tol = 1.e-200, Reduce = False):
    """Performs Gauss-Jordan elimination on matrix M to put 
it in reduced row echelon form.
Assumes last column is a column of constants."""
    lead = 0
    rCount, cCount = M.shape
    cCount = cCount - 1
    for r in range(rCount):
        if lead >= cCount:
            if Reduce:
                M = ReduceMatrix(M, Tol)
            return M
        i = r
        while abs(M[i,lead]) <= Tol:
            i += 1
            if i == rCount:
                i = r
                lead += 1
                if lead == cCount:
                    if Reduce:
                        M = ReduceMatrix(M, Tol)
                    return M
        M[i,:], M[r,:] = M[r,:].copy(), M[i,:].copy()
        M[r,:] = M[r,:] / M[r,lead]
        for i in range(rCount):
            if i != r:
                M[i,:] = M[i,:] - M[i,lead] * M[r,:]
        lead += 1
    if Reduce:
        M = ReduceMatrix(M, Tol)
    return M
    

def DependentVars(M, Tol = 1.e-200, IndepPriorities = []):
    """Returns a list of (var index, coefs, const) for expressing
dependent variables as linear functions of other variables, given the system
of linear equations described by the matrix M.  IndepPriorities
is a list of variable indices desired to be independent."""
    #reorder M by the priorities, putting desired independent variables at the end
    NewM = M.copy()
    CopiedCol = []
    CurNewCol = M.shape[1] - 2
    Map1 = np.zeros(M.shape[1] - 1, dtype=int)
    Map2 = np.zeros(M.shape[1] - 1, dtype=int)
    #start copying columns to the new matrix
    for i in IndepPriorities + range(M.shape[1] - 2, -1, -1):
        #skip if we already did this column
        if i in CopiedCol: 
            continue
        #sanity check on i
        if i < 0 or i > M.shape[1]:
            raise ValueError("Invalid variable index %d for matrix M." % i)
        #copy to last column of new matrix
        NewM[:,CurNewCol] = M[:,i]
        #make a map for converting to new column layout; value is index of new
        Map1[i] = CurNewCol
        #make a map for converting to old column layout; value is index of old
        Map2[CurNewCol] = i
        #record copy and update counters
        CopiedCol.append(i)
        CurNewCol -= 1
        if CurNewCol < 0:
            break
    #now perform GJ elimination
    NewM = GaussJordan(NewM, Tol = Tol, Reduce = True)
    #loop through linear constraints to create return array
    ret = []
    for Row in NewM:
        for (i, v) in enumerate(Row[:-1]):
            if abs(v) <= Tol:
                continue
            else:
                Ind = Map2[i]
                Coef = -Row[:-1].copy()
                Coef[i] = 0.
                Coef = Coef.take(Map1)
                Const = Row[-1]
                ret.append((Ind, Coef, Const))
                break
    #return the results
    return ret


