"""Functions to construct the EDS grid from simulated data"""

import numpy as np

def Dist( X,Xi ):
    """Distance from Xi to all points in X
    X is a m x n array
    Xi is a 1-d array of length n"""

    return np.sqrt(np.sum( (X - Xi)**2,axis = 1))


def MEpsGrid(P,M, tol,eps0,order=[],FirstRun = True):
    """Looks for eps such that EpsGrid returns M points

    (grid, eps) = MEpsGrid(P,M, tol,eps0,order=[],M0=[])

    Inputs: P is n x d array of simulated data
    M is desired number of points in the grid
    tol is tolerance around M
    eps0  2-tuple that brackets the distance between points
    order is the order of the points to eliminate.  Optional, default generates it at random.
    FirstRun (optional defaults true) indicates that we should check the boundaries.
    Outputs: (grid, eps) where grid is a M' x d grid  and eps is the distance between points."""



    d = P.shape[1]

    if len(order) == 0:
        n = P.shape[0]
        np.random.seed(41311)
        order = np.argsort(np.random.rand(n))


    if FirstRun:
        # print "Constructing Eps-Grid"
        Pe = EpsGrid(P,eps0[0],order)
        a = Pe.shape[0]
        # print "a = %s" % a
        if a >= M-tol and a <= M+tol:
            return (Pe, eps0[0])
        if a < M:
            raise RuntimeError("eps0[0] too large")



        Pe = EpsGrid(P,eps0[1],order)
        b = Pe.shape[0]
        # print "b = %s" % b
        if b >= M-tol and b <= M+tol:
            return (Pe, eps0[1])
        if b > M:
            raise RuntimeError("eps0[1] too small")



    eps = np.mean(eps0)
    Pe = EpsGrid(P,eps,order)
    c = Pe.shape[0]
    # print "c = %s" % c

    if c >= M-tol and c <= M+tol:
        return (Pe,eps)

    if c > M:  #midpoint produces too many points so it is too small
        epsnew = (eps, eps0[1])
        (Pe,eps) = MEpsGrid(P,M, tol,epsnew,order,False)
    else: # midpoint produces too few points so it is too large
        epsnew = (eps0[0], eps)
        (Pe,eps) = MEpsGrid(P,M, tol,epsnew,order,False)


    return (Pe, eps)




def EpsGrid(P,eps,order):
    """Create an Epsilon grid for data in P.
    Pe = EpsGrid(P,eps,order)
    Inputs: P is a n x d array where each row is a point in the state space.
    eps is the scalar distance such that there are no points in the grid closer than eps from other points.
    order is length-n array that gives the order in which the points in P are included in the grid.
    Outputs: Pe is the m x d array of grid points."""

    (n, d) = P.shape
    Pe = np.zeros((n,d))

    NotDeleted = np.ones(n, dtype=bool)

    ctr = 0
    for i in range(n):
        if NotDeleted[order[i]]:
            Pe[ctr,:] = P[order[i],:]
            ctr +=  1
            dstnce = Dist(P[NotDeleted,:],P[order[i],:])
            NotDeleted[NotDeleted] = dstnce > eps


        if not NotDeleted.any():
            break


    Pe = Pe[0:ctr,:]


    return Pe



def GetGrid(state,thin,M,tol,ElimLowProbPoints = True):
    """Takes simulated data and generates an EDS grid.

    (grid, eps, PCgrid) = GetGrid(state,thin,M,tol)

    Inputs: state is a d x T array of simulate states.
    thin is the fraction to keep after thinning
    M is scalar target number of grid points

    Outputs:
    grid is the M+/-tol x d grid
    eps is the distance between grid points
    PCgrid is the grid in principal component space"""


    (nState, simT) = state.shape

    #thin to [thin]%
    EDSSIMt = int(thin*simT)
    # print "EDSSIMt = %s" % EDSSIMt

    state = state[:,::int(1.0/thin)].transpose()

    # pricnipal components
    statemean = np.mean(state,axis = 0)
    ss = state - statemean
    statestd = np.std(state,axis = 0)
    if any(statestd == 0):
        raise RuntimeError('Zero standard deviation encountered')
    ss = np.dot(ss,np.diag(1./statestd))


    U, S, V = np.linalg.svd(ss)
    V = V.transpose()

    PC = np.dot(ss,V)
    D = np.std(PC,axis=0)
    PC = np.dot(PC , np.diag(1.0/D))
    PC2ss = np.dot(np.diag(D),V.transpose())

    tmp = np.corrcoef(PC,rowvar = 0) - np.eye(PC.shape[1])
    #assert(abs(tmp).all() < 1e-3)

    if ElimLowProbPoints:
        # eliminate low probability points
        # print "EDS -- computing densities"
        EDSh = EDSSIMt**(-1.0/(nState+4))
        EDSdelta = 0.1  # fraction to drop

        dens = np.zeros(EDSSIMt)

        # print "PC.shape, dens.shape"
        # print PC.shape
        # print dens.shape

        for i in range(EDSSIMt):
            dens[i] = (1.0/(EDSSIMt*(2*np.pi)**(nState/2.0)*EDSh**nState))*np.sum(np.exp(-Dist(PC,PC[i,:]))/(2*EDSh**2))

        I = np.argsort(dens)

        #discard lowest 10%
        ndisc = int(np.round(EDSSIMt*EDSdelta))
        PC = PC[I[ndisc+1:],:]

    # EDS 4 - create eps grid
    eps0 = (0.01, 10.)
    (PCgrid, eps) = MEpsGrid(PC,M, tol, eps0)


    Mp = PCgrid.shape[0]
    grid = np.dot(PCgrid,np.dot(PC2ss,np.diag(statestd))) + statemean


    return (grid, eps, PCgrid)
