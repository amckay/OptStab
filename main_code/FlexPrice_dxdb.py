"""
Unpacks the numerical optimization over b and tau using the propositions
making use of the no-mortality simplification.

Results are saved, plot them with Unpack_Plot.ipynb

Usage:
python OptStab_Unpack.py -param {b,tau} -mode {fig,point}

"""
import numpy as np
from scipy.linalg import block_diag
import os
import copy
import sys
from multiprocessing import Pool
from functools import partial
import ModelFramework.PolyInt as PI
from ModelFramework import EDS
from ModelFramework.Drivers import SobolGrid
from ModelFramework.Support import NaNInfError
import time

from OptStab_FlexiblePrices import (bmin, bmax, taumin, taumax, taustep, OptStabDynamics,IncProc,b,tau,
                    VarNums,shock_covar,nquad,WelfareStates,DynamicsStates,poly_ord,
                    adexp, adlog, exp, log, beta, upsilon,psi_1,psi_2,kappa,psy,gamma,gamma_0,delta,chi,
                    ieG,iA, iMU, iY, ih, iM, iu, iulag, iq,iEAlph,iElogAlpha,iVn,iJ,
                    T,nGrid,VERSION, bstep, VarNums)

bcalib = b
taucalib = tau



# globals
OutDir = os.path.join(os.path.pardir,os.path.pardir,'data','results')
NWORKERS = 4
NCOMP = 7

#order of Envelope conditions ulag, EAlph, ElogAlpha, ulag_given_x
nEnv = 4
iEnv_ulag, iEnv_EAlph, iEnv_ElogAlpha, iEnv_ulag_given_x = range(nEnv)



with np.load(os.path.join(OutDir,'SearchResults_v8.npz')) as X:
    bOpt, tauOpt = X['OptCyc']


NREP = 1000
Tsim = 600




def SolveMod(b,tau):
    # Create the welfare model object
    M2 = OptStabDynamics(VarNums,shock_covar,PI.Poly(poly_ord),IncProc,b,tau,nquad,WelfareStates)

    M2.exp = adexp
    M2.log = adlog
    M2.Perturb(SetF = True)
    M2.exp = exp
    M2.log = log
    # M2.SolveEDS(T,nGrid,damp = [0.995,0.995],verbose = False, PlotGrids = False)
    return M2




def Getdxdparam(Mda,Mdb,state):
    """Returns a finite difference derivative of x = M with respect to param"""


    Xa = np.zeros((Mda.nXYZ,state.shape[1]))

    Xa[Mda.nX:Mda.nXY] = Mda.F(Xa[Mda.interpstates])
    Xa[Mda.nXY:] = Mda.Static(Xa)

    Xb = Xa.copy()
    #Xb[iulag] = Xa[iulag] + (1-Xa[iq]*Xa[iM])*(Mdb.ubar-Mda.ubar)
    Xb[Mdb.nX:Mdb.nXY] = Mdb.F(Xb[Mdb.interpstates])
    Xb[Mdb.nXY:] = Mdb.Static(Xb)


    D = Mdb.b - Mda.b

    return (Xb[iM] - Xa[iM])/D



if __name__ == '__main__':



    # create a grid
    interpstates = WelfareStates
    F = PI.Poly(poly_ord)
    M = OptStabDynamics(VarNums,shock_covar,F,IncProc,(bmin+bmax)/2,(taumin+taumax)/2,nquad,interpstates)
    M.TestSteadyState()

    # --- first order perturbation --
    M.exp = adexp
    M.log = adlog
    PertSol = M.Perturb(SetF = True)
    M.exp = exp
    M.log = log

    # -- create a grid --
    Xsim = M.EDSSimulate(T)
    noninterpstates = np.setdiff1d(range(M.nX),M.interpstates)
    THIN = max(0.1, 2 * float(nGrid)/float(T))
    EDS_points_tol = 2

    grid_a, eps, PCGrid = EDS.GetGrid(Xsim[M.interpstates],THIN,nGrid,EDS_points_tol)
    grid_a = grid_a.T
    grid = np.zeros((M.nX,grid_a.shape[1]))
    grid[M.interpstates] = grid_a
    grid[noninterpstates] = Xsim[noninterpstates,:grid.shape[1]]




    paramgrid = np.arange(bmin,bmax+bstep,bstep)
    MList = map(lambda bb: SolveMod(bb,tauOpt),paramgrid)

    tst = np.array([z.F.coef[0,:2] for z in MList])
    tst = np.dot(tst,grid[:2,:])
    dxdb = np.diff(tst,axis = 0) / np.diff(paramgrid)[np.newaxis].T
    print("Max |dxdb| = {}".format(np.max(np.abs(dxdb))))
