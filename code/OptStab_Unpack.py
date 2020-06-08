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

OutDir = os.path.join(os.path.pardir,'data','results')

from OptStab import (bmin, bmax, taumin, taumax, taustep, OptStabDynamics,IncProc,b,tau,
                    VarNums,shock_covar,nquad,WelfareStates,poly_ord,
                    adexp, adlog, exp, log, beta, upsilon,psi_1,psi_2,kappa,psy,gamma,gamma_0,delta,chi,SSGRatio,zeta_2,
                    ieG,iA, iMP, iMU, iY, iSA, ih, iM, iu, iulag, iq,iEAlph,iElogAlpha,idVndb,idVndtau,iVn,iJ,iS,
                    T,nGrid,VERSION,CLArgs)


print("\nOptions:")
print("mode = {}".format(CLArgs.mode))
print("param = {}".format(CLArgs.param))
print("")

bcalib = b
taucalib = tau



# globals
NWORKERS = 4
NCOMP = 7

#order of Envelope conditions ulag, EAlph, ElogAlpha, S, ulag_given_x
nEnv = 5
iEnv_ulag, iEnv_EAlph, iEnv_ElogAlpha, iEnv_S, iEnv_ulag_given_x = range(nEnv)


if CLArgs.mode == "fig":
    OptVer = VERSION
else:
    OptVer = 1

with np.load(os.path.join(OutDir,'SearchResults_v'+str(OptVer)+'.npz')) as X:
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



def SkillDeriv(Md,X):
    """Gives d/du of  E_i eps^(1-tau)  and E_i log( eps) """

    if VERSION == 16:
        return np.zeros(Mom1.shape)
    elif VERSION == 31:
        utemp = 1.25*(X[iu]-Md.ubar)+1.0*Md.ubar
    else:
        utemp = X[iu]

    Mom1 = Md.IP.get_Moments(utemp,Md.ubar,Md.tau)
    Mom1 = np.array((Mom1[0], Mom1[1]  ))
    eps = 0.001
    Mom2 = Md.IP.get_Moments(utemp+eps,Md.ubar,Md.tau)
    Mom2 = np.array((Mom2[0], Mom2[1]  ))
    return (Mom2-Mom1)/eps


if CLArgs.param == "b":

    def Insurance(Md,X):
        """Gives insurance component u/b - u/(1-u+ub)"""
        u = X[iu]
        b = Md.b()
        return u/b - u/(1-u+u*b)

    def Incentives(Md,X,EEnvPrime):
        """Gives incentives component (dW/du)(du/db|_M)-upsilon q**kappa  (dq/db|_M)"""
        u = X[iu]
        b = Md.b()



        CplusG = 1./X[iMU] * (1+SSGRatio*X[ieG])

        SD = SkillDeriv(Md,X)
        SkillRisk =  beta * (  EEnvPrime[iEnv_EAlph] *(1-delta)*X[iEAlph]*SD[0]
                            +  EEnvPrime[iEnv_ElogAlpha]*(1-delta)*SD[1])

        dWdu = (log(b) + (1-b)/(1-u+u*b) - (1+chi)*X[iA]*X[ih]/CplusG/X[iSA]
                + (1+chi)*psi_1*X[iM]**psi_2/CplusG+gamma_0*X[ih]**(1+gamma)/(1+gamma)-psy + SkillRisk)
        dqdb_M =  (1./kappa)*X[iq] * X[idVndb] / X[iVn]
        dudq_M = -X[iM]*(X[iulag]+upsilon*(1-X[iulag]))
        dudb_M = dudq_M * dqdb_M

        dWdh =  ( (1+chi)*X[iA]/CplusG/X[iSA] - gamma_0 * X[ih]**gamma ) * (1-u)

        # dhdq = -zeta_2/(1+gamma)*X[ih]/(1-X[iu])*dudq_M
        # dhdubar = zeta_2/(1+gamma)*X[ih]/(1-Md.ubar)
        # dhdb =   dhdq/kappa*X[iq]/X[iVn]*X[idVndb] + dhdubar * Md.dubardb

        dhdu = -zeta_2/(1+gamma)*X[ih]/(1-X[iu])


        XSS = Md.XSS
        dhdubar = zeta_2/(1+gamma)*X[ih]/(1-XSS[iu])
        dubardq_M = -XSS[iM]*(XSS[iulag]+upsilon*(1-XSS[iulag]))
        dqbardb_M = (1./kappa)*XSS[iq] * XSS[idVndb] / XSS[iVn]


        dhdb =   dhdu*dudq_M*dqdb_M + dhdubar * dubardq_M*dqbardb_M


        return dWdu * dudb_M - (X[iulag]+upsilon*(1-X[iulag]))*X[iM]*X[iVn]*dqdb_M  + dWdh * dhdb + beta*dudb_M*EEnvPrime[iEnv_ulag]

else:
    def Insurance(Md,X):
        """Gives insurance component.  This is only the part related to reducing the impact of future shocks
        and not the benefit of redistributing across already existing dispersion in alpha."""
        if VERSION == 16:
            utemp = 0.0*X[iu]+1.0*Md.ubar
        elif VERSION == 31:
            utemp = 1.25*(X[iu]-Md.ubar)+1.0*Md.ubar
        else:
            utemp = X[iu]

        Mom = Md.IP.get_Moments(utemp,Md.ubar,Md.tau)
        return beta*(-Mom[1]+Mom[3]/Mom[0])

    def Incentives(Md,X,EEnvPrime):
        """Gives incentives component  """

        u = X[iu]
        b = Md.b()
        CplusG = 1./X[iMU] * (1+SSGRatio*X[ieG])

        SD = SkillDeriv(Md,X)
        SkillRisk =  beta * (  EEnvPrime[iEnv_EAlph] *(1-delta)*X[iEAlph]*SD[0]
                            +  EEnvPrime[iEnv_ElogAlpha]*(1-delta)*SD[1])

        dWdu = (log(b) + (1-b)/(1-u+u*b) - (1+chi)*X[iA]*X[ih]/CplusG/X[iSA]
                + (1+chi)*psi_1*X[iM]**psi_2/CplusG+gamma_0*X[ih]**(1+gamma)/(1+gamma)-psy  + SkillRisk)
        dWdh =  ( (1+chi)*X[iA]/CplusG/X[iSA] - gamma_0 * X[ih]**gamma ) * (1-u)

        dudq = -X[iM]*(X[iulag]+upsilon*(1-X[iulag]))
        dhdq = -zeta_2/(1+gamma)*X[ih]/(1-X[iu])*dudq
        dhdubar = zeta_2/(1+gamma)*X[ih]/(1-Md.ubar)

        XSS = Md.XSS
        dhdubar = zeta_2/(1+gamma)*X[ih]/(1-XSS[iu])
        dubardq_M = -XSS[iM]*(XSS[iulag]+upsilon*(1-XSS[iulag]))
        dqbardtau_M = (1./kappa)*XSS[iq] * XSS[idVndtau] / XSS[iVn]

        dhdtau = - X[ih]/(1+gamma)/(1-Md.tau) + dhdq/kappa*X[iq]/X[iVn]*X[idVndtau] + dhdubar * dubardq_M*dqbardtau_M

        dqdtau = (1./kappa)*X[iq] * X[idVndtau] / X[iVn]
        dudtau = dudq * dqdtau

        return dWdh*dhdtau + dWdu * dudtau  - (X[iulag]+upsilon*(1-X[iulag]))*X[iM]*X[iVn]*dqdtau + beta*dudtau*EEnvPrime[iEnv_ulag]



    def GetRedistributionValue(Md):
        """Only relevnat for -param tau
        Marginal insurance/redistribution value of increasing tau from smoothing
        already existing dispersion in alpha"""
        X = Md.SteadyState()
        Mom = Md.IP.get_Moments(Md.ubar,Md.ubar,Md.tau)
        Ealpha1taulogalpha = (1-delta)*X[iEAlph]*Mom[3]/(1-(1-delta)*Mom[0])
        return (-X[iElogAlpha]+Ealpha1taulogalpha/X[iEAlph] )/(1-beta)


def MacroStabilization(Md,X,EEnvPrime):
    """Gives an array where cols are components of dW/dx
    Cols correspond to: labor wedge, price dispersion, Hosios, unemployment risk, skill risk"""

    dXdx = GetdXdx(Md,X)
    u = X[iu]
    b = Md.b()
    CplusG = 1./X[iMU] * (1+SSGRatio*X[ieG])


    Hosios = ((1+chi)/CplusG* (-X[iA]/X[iSA]*X[ih] + psi_1*X[iM]**psi_2) * dXdx[iu]
                - (1+chi)/CplusG * psi_1*psi_2*X[iM]**(psi_2-1)*(1-X[iu]-(1-upsilon)*(1-X[iulag]))
                + beta * dXdx[iu] * EEnvPrime[iEnv_ulag] )

    LaborWedge = (1-u)*((1+chi)*X[iA]/CplusG/X[iSA] - X[ih]**gamma ) * dXdx[ih]

    PriceDispersion = -(1+chi)*X[iY]/X[iSA]/CplusG   * dXdx[iSA] + beta * dXdx[iSA] * EEnvPrime[iEnv_S]

    UnemploymentRisk = ( (log(b)+gamma_0*X[ih]**(1+gamma)/(1+gamma) - psy + (1-b)/(1-u+u*b) ) * dXdx[iu]
                        -(X[iulag] + upsilon*(1-X[iulag])) * X[iq]**kappa * dXdx[iq]  )


    SD = SkillDeriv(Md,X)*dXdx[iu]
    SkillRisk =  beta * (  EEnvPrime[iEnv_EAlph] *(1-delta)*X[iEAlph]*SD[0]
                        +  EEnvPrime[iEnv_ElogAlpha]*(1-delta)*SD[1])


    return np.hstack((LaborWedge,PriceDispersion,Hosios,UnemploymentRisk,SkillRisk))


def GetdXdx(Md,X):
    """Gives the derivative of X with respect to x = M.
    To construct this derivative, we use a monetary shock."""
    assert X.shape[1] == 1
    Xa = X.copy()
    Xb = X.copy()

    eps = -0.005
    Xb[iMP] = Xb[iMP] + eps


    Xa[Md.nX:Md.nXY] = Md.F(Xa[Md.interpstates])
    Xa[Md.nXY:] = Md.Static(Xa)

    Xb[Md.nX:Md.nXY] = Md.F(Xb[Md.interpstates])
    Xb[Md.nXY:] = Md.Static(Xb)

    return (Xb-Xa)/(Xb[iM] - Xa[iM])


def Getdxdparam(Mda,Mdb,Xa):
    """Returns a finite difference derivative of x = M with respect to param"""

    Xb = Xa.copy()
    #Xb[iulag] = Xa[iulag] + (1-Xa[iq]*Xa[iM])*(Mdb.ubar-Mda.ubar)
    Xb[Mdb.nX:Mdb.nXY] = Mdb.F(Xb[Mdb.interpstates])
    Xb[Mdb.nXY:] = Mdb.Static(Xb)

    if CLArgs.param == "b":
        D = Mdb.b() - Mda.b()
    else:
        D = Mdb.tau - Mda.tau

    return (Xb[iM] - Xa[iM])/D



def Dynamics(Md,X,shockprime):
    """Returns the next value of X given next periods shocks shockprime"""
    assert X.shape[1] == 1

    Xp = X.copy()
    Xp[:Md.nX] = Md.StateTrans(X,shockprime)
    Xp[Md.nX:Md.nXY] = Md.F(Xp[Md.interpstates])
    Xp[Md.nXY:] = Md.Static(Xp)

    return Xp

def GetdxdState(Md,X,iState):
    #compute derivatives wrt state number iState
    Y = X.copy()
    eps = 0.0005
    Y[iState] += eps
    Y[Md.nX:Md.nXY] = Md.F(Y[Md.interpstates])
    Y[Md.nXY:] = Md.Static(Y)
    return (Y-X)/eps

def rshp(Z,nBP):
    return np.tile( Z[np.newaxis].T, (1,nBP))

def GetPD(Md,X,nBP,iState):
    '''Helper function for GetExpectedEnvCond.  Returns period contribution P and Discounting matrix D'''

    dX = GetdxdState(Md,X,iState)

    u = X[iu]
    b = Md.b()
    CplusG = 1./X[iMU] * (1+SSGRatio*X[ieG])

    if VERSION == 16:
        utemp = 0.0*X[iu]+1.0*Md.ubar
    elif VERSION == 31:
        utemp = 1.25*(X[iu]-Md.ubar)+1.0*Md.ubar
    else:
        utemp = X[iu]

    eps1mt = Md.IP.get_Moments(utemp,Md.ubar,Md.tau)[0]


    dWdu = (log(b) + (1-b)/(1-u+u*b) - (1+chi)*X[iA]*X[ih]/CplusG/X[iSA]
        +gamma_0*X[ih]**(1+gamma)/(1+gamma)-psy )
    dWdJ = (1+chi) /CplusG
    dWdh = ((1+chi)/CplusG*X[iA]/X[iSA] -gamma_0*X[ih]**gamma)*(1-u)
    dWdq = (X[iulag] + upsilon*(1-X[iulag])) * X[iq]**kappa
    dWdS = -(1+chi)/CplusG*X[iY]/X[iSA]

    deps1mtdu, dlogepsdu = SkillDeriv(Md,X)


    P = (     dWdu * dX[iu]
            + dWdJ * dX[iJ]
            + dWdh * dX[ih]
            + dWdq * dX[iq]
            + dWdS * dX[iSA]   )

    D = beta * np.hstack((  rshp( dX[iu] ,nBP),
                     rshp( (1-delta)*eps1mt*dX[iEAlph] + (1-delta)*X[iEAlph]*deps1mtdu*dX[iu] ,nBP),
                     rshp( (1-delta)*dX[iElogAlpha] + (1-delta)*dlogepsdu*dX[iu] ,nBP),
                     rshp( dX[iSA] ,nBP),
                     np.zeros((X.shape[1],nBP))  ))

    return P, D


def GetExpectedEnvCond(Md,grid):
    #order ulag, EAlph, ElogAlpha, S, ulag_given_x


    X = np.zeros((Md.nXYZ,grid.shape[1]))
    X[:Md.nX] = grid
    X[Md.nX:Md.nXY] = Md.F(X[Md.interpstates])
    X[Md.nXY:] = Md.Static(X)

    BP = PI.PolyIntBasis(grid,Md.F.ord).T
    nBP = BP.shape[1]
    EBPPrime = 0.0
    for quad_i in range(np.prod(Md.quad['n'])):
        Stateprime = Md.StateTrans(X,Md.quad['x'][quad_i])
        EBPPrime += Md.quad['w'][quad_i] * PI.PolyIntBasis(Stateprime,Md.F.ord).T



    m = grid.shape[1]
    n = m * nEnv
    P = np.zeros(n)
    D = np.zeros((n,nEnv*nBP))


    P[:m], D[:m] =  GetPD(Md,X,nBP,iulag)
    P[m:2*m], D[m:2*m] =  GetPD(Md,X,nBP,iEAlph)
    P[2*m:3*m], D[2*m:3*m] =  GetPD(Md,X,nBP,iElogAlpha)
    P[3*m:4*m], D[3*m:4*m] =  GetPD(Md,X,nBP,iS)


    P[:m] += -(1-upsilon)*X[iq]**(1+kappa)/(1+kappa)
    P[m:2*m] += - 1./X[iEAlph]
    P[2*m:3*m] += 1. - Md.tau
    P[3*m:4*m] +=  0.

    # the last one is envelope wrt ulag holding x fixed
    CplusG = 1./X[iMU] * (1+SSGRatio*X[ieG])
    P[4*m:] = ( (1+chi)/CplusG *( ( -X[iA] * X[ih] / X[iSA] + psi_1 * X[iM]**psi_2) * (1-upsilon)*(1-X[iq]*X[iM])  - psi_1 * X[iM]**psi_2*(1-upsilon))
                + (log(Md.b())+(1-Md.b())/(1-X[iu]*(1-Md.b())) +  gamma_0 * (X[ih]**(1+gamma)/(1+gamma)-psy)) * (1-upsilon)*(1-X[iq]*X[iM]) - (1-upsilon)*X[iq]**(1+kappa)/(1+kappa) )

    D[4*m:,4*nBP:] = beta * rshp( (1-upsilon)*(1-X[iq]*X[iM]) ,nBP)

    coefs = np.linalg.lstsq(block_diag( *(nEnv*(BP,)) ) - D * np.tile(EBPPrime,(nEnv,nEnv)),P,rcond=None)[0]
    coefs = coefs.reshape(nEnv,nBP).T
    ExpectedEnvCondOnGrid = np.dot(EBPPrime, coefs).T

    FLEnv = PI.Poly(2)
    FLEnv.FitXY(grid , ExpectedEnvCondOnGrid)
    return FLEnv


def ExpectedMarginalWelfare(M,grid,i):
    np.random.seed(48230)
    y = np.zeros((NREP,NCOMP+1,Tsim))
    mu = np.zeros((NREP,Tsim,2))
    EEnvPrime = np.zeros((nEnv,1))

    # -- fit polynomial to expectations --
    # interpolate the forward looking expectation of envelope conditions
    if i < len(M)-1:
        FLEnv = GetExpectedEnvCond(M[i],grid)

    for it in range(NREP):

            nshock = M[i].shock_covar.shape[0]
            eps = np.dot(np.linalg.cholesky(M[i].shock_covar),np.random.randn(nshock,Tsim))

            X = M[i].SteadyState()[np.newaxis].T
            try:
                for t in range(Tsim):
                    if i < len(M)-1:
                        dxdparam = Getdxdparam(M[i],M[i+1],X)
                        EEnvPrime[:] = FLEnv(X[:M[i].nX])
                        Insur = Insurance(M[i],X)
                        Incen = Incentives(M[i],X,EEnvPrime)
                        MS = MacroStabilization(M[i],X,EEnvPrime)
                        y[it,:-1,t] = np.hstack((Insur,Incen,MS))
                        y[it,-1,t] = dxdparam
                    mu[it,t,:] = X[[iM,iu]].flatten()
                    X = Dynamics(M[i],X,eps[:,t])

            except NaNInfError:
                print("Bad values in iteration {0} for b = {1}, tau = {2}".format(it,M[i].b(),M[i].tau))



    if i < len(M)-1:
        if CLArgs.param == "tau":
            # add redistribution to insurance benefit
            y[:,0,0] += GetRedistributionValue(M[i])

        return y, mu
    else:
        return mu


if __name__ == '__main__':
#    from numba import jit

#    @jit
    def PostProcess(EMWRaw,muList):
        """Take expectations and sums across t"""
        KeepIt = np.ones(NREP,dtype=bool)
        for res_b in EMWRaw:
            for (i,res_rep) in enumerate(res_b):
                if np.isnan(res_rep).any():
                    KeepIt[i] = False

        dxdparam_2 = (muList[1:,:,:,0] - muList[:-1,:,:,0]) / np.diff(paramgrid).reshape(len(paramgrid)-1,1,1)


        EMW             = np.zeros((len(paramgrid)-1,NCOMP))
        EdWdxEdxdparam  = np.zeros((len(paramgrid)-1,NCOMP))
        Cov             = np.zeros((len(paramgrid)-1,NCOMP))
        Cov_2             = np.zeros((len(paramgrid)-1))
        std_dWdx        = np.zeros((len(paramgrid)-1))
        std_dxdparam    = np.zeros((len(paramgrid)-1))
        for (i,res_b) in enumerate(EMWRaw):
            u = muList[i][KeepIt,:,1]
            IIdWdx = res_b[KeepIt,:-1,:]
            dxdparam = res_b[KeepIt,-1,:]
            dxdparam_2i = dxdparam_2[i,KeepIt,:]
            tmp = np.ones((dxdparam.shape[0],NCOMP,dxdparam.shape[1])) # don't multiply the insurance and incentives terms against dxdparam
            tmp[:,2:,:] = np.tile(np.reshape(dxdparam,(dxdparam.shape[0],1,dxdparam.shape[1])),(1,NCOMP-2,1))
            EMW[i,:] = (IIdWdx * tmp * beta**np.arange(Tsim)).sum(axis=2).sum(axis=0)/ KeepIt.sum()
            EdWdxEdxdparam[i,:] = (IIdWdx.sum(axis=0)/ KeepIt.sum() * tmp.sum(axis=0)/ KeepIt.sum() * beta**np.arange(Tsim)).sum(axis=1)
            for t in range(Tsim):
                std_dWdx[i] += beta**t * np.std(IIdWdx[:,2:,t].sum(axis=1))
                std_dxdparam[i] += beta**t * np.std(dxdparam[:,t])
                for j in range(2,NCOMP):
                    Cov[i,j] += beta**t * np.cov(IIdWdx[:,j,t].flatten(),dxdparam[:,t].flatten())[0,1]
                Cov_2 += beta**t * np.cov(IIdWdx[:,2:,t].sum(axis=1), dxdparam_2i[:,t].flatten() )[0,1]

        return EMW, EdWdxEdxdparam, Cov, std_dWdx, std_dxdparam, Cov_2




    # load the saved grid
    fname = os.path.join(OutDir,'GridBounds_v'+str(VERSION)+'.npz')
    with np.load(fname) as X:
        savedgrid = X['savedgrid']

    grid_a = EDS.GetGrid(savedgrid[WelfareStates],1.0,nGrid,2,ElimLowProbPoints = False)[0].T
    grid = np.zeros((VarNums[0],grid_a.shape[1]))
    grid[WelfareStates] = grid_a
    noninterpstates = np.setdiff1d(range(VarNums[0]),WelfareStates)
    grid[noninterpstates] = savedgrid[noninterpstates,:grid.shape[1]]

    # add some Sobol points
    GridWidth = 1.5 * ( grid.max(axis = 1) - grid.min(axis = 1) )
    Sob = SobolGrid(grid.mean(axis = 1),GridWidth, nSobol = 20)
    grid = np.hstack((grid,Sob))



    # # -- main work --
    # Solve the model
    if CLArgs.mode == "fig":
        if CLArgs.param == "b":
            bstep = 0.016
            paramgrid = np.arange(bmin,bmax+bstep,bstep)
            MList = map(lambda bb: SolveMod(bb,tauOpt),paramgrid)
        else:
            paramgrid = np.arange(taumin,taumax+taustep,taustep)
            MList = map(lambda t: SolveMod(bOpt,t),paramgrid)



    else:  # point
        if CLArgs.param == "b":
            paramgrid = np.array([bOpt,bOpt+0.01])
            MList = map(lambda bb: SolveMod(bb,tauOpt),paramgrid)
        else:
            paramgrid = np.array([tauOpt,tauOpt+0.01])
            MList = map(lambda t: SolveMod(bOpt,t),paramgrid)





    # do the unpacking
    start_time = time.time()
    if CLArgs.mode == "fig":
        p = Pool(NWORKERS)
        mymap = p.map
    else:
        mymap = map

    EMWRaw_u = mymap( partial(ExpectedMarginalWelfare,MList,grid),  range(len(paramgrid)-1) )
    EMWRaw = np.array([x[0] for x in EMWRaw_u])
    mu = [x[1] for x in EMWRaw_u]
    mu.append(ExpectedMarginalWelfare(MList,grid,len(paramgrid)-1))
    mu = np.array(mu)


    print 'Elapsed time = {}'.format(time.time()-start_time)
    # np.save('/projectnb/autostab/Data/scratch/EMWRaw.npz',EMWRaw)
    # np.save('/projectnb/autostab/Data/scratch/mu.npz',mu)
    start_time = time.time()
    EMW, EdWdxEdxdparam, Cov, std_dWdx, std_dxdparam, Cov_2 = PostProcess(EMWRaw, mu)
    print 'Elapsed time = {}'.format(time.time()-start_time)


    np.savez(os.path.join(OutDir,'Unpack_v'+str(VERSION)+'_param_'+CLArgs.param+'_mode_'+CLArgs.mode + '.npz'),
                ExpectedMarginalWelfare = EMW,
                EdWdxEdxdparam = EdWdxEdxdparam,
                Cov = Cov,
                Cov_2 = Cov_2,
                std_dWdx = std_dWdx,
                std_dxdparam = std_dxdparam,
                paramgrid = paramgrid)
