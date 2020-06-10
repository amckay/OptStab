import numpy as np
from ModelFramework.CompEcon import qnwnorm
import ModelFramework.Drivers as Drivers
import ModelFramework.EDS as EDS
import ModelFramework.PolyInt as PI
import matplotlib.pyplot as plt
import json
from numpy import exp, log
import ad
from scipy import optimize
import csv
from ModelFramework.AutoDiff import exp as adexp
from ModelFramework.AutoDiff import log as adlog
import time
import warnings
from OptStab_IncomeProcess import IncomeProcess
from scipy import optimize
import os

VERSION = 8

OutDir = os.path.join(os.path.pardir,os.path.pardir,'data','results')

# indexing
(iA, ieG, iEAlph, iElogAlpha, iulag,
 iM, iVn, iV_a, iV_b, iV_c, iV, iYhat,
 iG, iY, ih, iw, iu, iJ, iq, iMU) = range(20)

iVRange = np.array([ iV_a, iV_b, iV_c, iV])

execfile('OptStab_calibration.py')

class OptStabDynamics(Drivers.Model):

    def __init__(self,VarNums,shock_covar,F,IncProc,b,tau,nquad,interpstates):
        self.log = log
        self.exp = exp

        self.IP = IncProc
        self.Set_b_tau(b,tau)
        self.VarNums = VarNums


        # call the parent init
        super(OptStabDynamics,self).__init__(VarNums,shock_covar,F,nquad,interpstates)

    def Set_b_tau(self,b,tau):
        self.b = b
        self.tau = tau
        self.wbar = np.maximum(0.99*wcalib,wcalib * (b/bcalib)**WageElasticity)
        SS = self.SteadyState()
        self.hbar = SS[ih]
        self.ubar = SS[iu]
        self.Ybar = SS[iY]
        self.Jbar = SS[iJ]
        self.Mbar = SS[iM]


    def SteadyState(self):
            w = self.wbar
            A = 1.
            MP_shock = 1.
            G_shock = 1.
            etam = 1.0
            S = 1.0
            SA = S
            pi = 1.
            pstar = 1.


            # six equations in h, Y, H, q, J, u
            # guess Y, J, and Vn
            def g(X,Mode = 'Residuals'):
                Y, J, Vn = X  #unpack
                h = ( (1-self.tau) * w * S * Y/ ( A * (Y-J)  )) **(1./(1.+gamma))
                M = ((A/mu  - w)*h/psi_1/upsilon)**(1/psi_2)
                q = (M*Vn)**(1.0/kappa)
                u = upsilon*(1 - q * M) / (q*M + upsilon*(1 - q * M) )
                H = 1-u - (1-upsilon)*(1-u)
                # check Y, J, Vn
                Y2 = A * h * (1-u)
                J2 = psi_1 * M**psi_2 * H
                Vn2 = (-log(self.b) -   h**(1+gamma)/(1+gamma) + psy )/(1-beta*(1-upsilon)*(1-kappa/(1+kappa)*q*M) )
                if Mode.lower() == 'residuals':
                    return np.array((Y2 - Y, J2-J, Vn2-Vn))
                elif Mode.lower() == 'values':
                    return Y, J,  h, q, u, H, M, Vn
                else:
                    raise Exception('Unrecognized value for Mode')


            sol = optimize.root(g,(Ycalib,Jcalib,Vncalib),jac=False)
            Y, J, h, q, u, H, M, Vn = g(sol.x,Mode = 'values')

            IP_Moments = self.IP.get_Moments(u,u,self.tau)
            ElogAlpha = (1-delta)*IP_Moments[1]/delta
            EAlph = float(delta / (1.0 - (1.0 - delta)*IP_Moments[0]))

            EU = upsilon*(1-q*M)
            I = 1.0/ ( beta * ((1-EU) + EU/self.b) * IP_Moments[2])
            C = (Y-J)/(1+chi)
            G = chi * C
            pB = Y/(1-(1-theta)/I)
            pA = pB

            MU = 1.0/C

            dVndb = -1./self.b  / (1-beta*(1-upsilon)*(1-q*M))

            X = np.array((A, G_shock, EAlph, ElogAlpha, u,
                            M, Vn, 0.,0.,0.,0.,Y,
                            G, Y, h, w, u, J, q, MU ))


            V = np.dot(np.linalg.inv(np.eye(len(iVRange)) - np.diag(self.Discounting(X).flatten()) ) ,
                           self.PeriodUtility(X)    )

            X[iVRange] = V.flatten()

            return X




    def StateTrans(self,X,shockprime):
        """Order of states is:
           0    1      2       3          4
        ( A, G_shock, EAlph, ElogAlpha,  ulag,    )"""

        IP_Moments = self.IP.get_Moments(X[iu],self.ubar,self.tau)
        ElogEpsilonPrime = IP_Moments[1]

        return np.vstack((
            self.exp(rhoA * self.log(X[iA]) + shockprime[0] ),
            self.exp(rhoG * self.log(X[ieG]) + shockprime[1]),
            (1-delta) * X[iEAlph] * IP_Moments[0][0] + delta,
            (1-delta) * (X[iElogAlpha] + ElogEpsilonPrime),
            X[iu]
            ))

    def Static(self,X):
        """order of static variables is
          0   1   2    3   4    5  6    7
        (iG, iY, ih,  iw, iu, iJ, iq, iMU)"""

        if X[iVn].__class__ == np.ndarray and X[iVn][0].__class__ == np.float64  and any(np.isnan(X[iVn])):
            raise RuntimeError("X[iVn] NaN")

        if X[iM].__class__ == np.ndarray and X[iM][0].__class__ == np.float64  and any(np.isnan(X[iM])):
            raise RuntimeError("X[iM] NaN")

        # M  and Vn give q
        q = (X[iM]* X[iVn])**(1./kappa)
        if q.__class__ == np.ndarray and q[0].__class__ == np.float64  and any(np.isnan(q)):
            print any(X[iVn] < 0)
            print any(X[iM] < 0)
            raise RuntimeError("q NaN")
        # u transition, q, and M give u
        u = (1-q*X[iM])*(X[iulag]+upsilon*(1-X[iulag]))
        if u.__class__ == np.ndarray and u[0].__class__ == np.float64  and any(np.isnan(u)):
            raise RuntimeError("u NaN")
        # def J gives J
        J = psi_1 * X[iM]**psi_2 * ( (X[iulag]+upsilon*(1-X[iulag])) - u)
        # wage rule gives w
        w = self.wbar * X[iA] * ((1-u)/(1-self.ubar))**zeta_2
        if w.__class__ == np.ndarray and w[0].__class__ == np.float64  and any(np.isnan(w)):
            raise RuntimeError("w NaN")
        # labor supply gives h
        h = (((1-self.tau)*w) / (X[iA] * (1-J/X[iYhat])  ))**(1./(1.+gamma))
        if h.__class__ == np.ndarray and h[0].__class__ == np.float64  and any(np.isnan(h)):
            print any(w < 0.)
            raise RuntimeError("h NaN")

        # prod func gives Y
        Y = X[iA]*h*(1-u)
        # Resource constraint and G rule give C and G
        C = (Y-J)/(1+chi*X[ieG])
        G = chi * C * X[ieG]



        return np.vstack(( G, Y, h, w, u, J, q, 1./C))


    def PeriodUtility(self,X):
        """ Period contributions to welfare and envelope conditions for
        iVn, iV_a, iV_b, iV_c, iV, iEnv_EAlph, iEnv_ElogAlpha, iEnv_S, iEnv_ulag_Hos, iEnv_ulag_Risk"""

        if np.any((-log(self.b)- gamma_0 * (X[ih]**(1+gamma)/(1+gamma) -psy)).ravel() ) < 0:
            warnings.warn('Negative search effort encountered. A')

        if np.any(X[iq].ravel() ) < 0:
            warnings.warn('Negative search effort encountered. B')

        if X[iq].__class__ == np.ndarray and (X[iq][0].__class__ == np.float64 or X[iq][0].__class__ == np.ndarray ) and any(np.isnan(X[iq])):
            warnings.warn('NaN search effort encountered.')

        if np.any(X[ih].ravel() ) < 0:
            warnings.warn('Negative work effort encountered.')


        numsearchers = upsilon*(1-X[iulag]) + X[iulag]

        PU_VComponents = np.vstack( (  (1-self.tau)*X[iElogAlpha] - self.log(X[iEAlph]),
                                    X[iu] * self.log(self.b) - self.log(1- X[iu]*(1-self.b)),
                                    -self.log(np.maximum(X[iMU],0.001)) - (1-X[iu])*gamma_0 * (X[ih]**(1+gamma)/(1+gamma)-psy) - numsearchers * X[iq]**(1+kappa)/(1+kappa) + chi * self.log(np.maximum(X[iG],0.001))
                                    ))


        return np.vstack((PU_VComponents, np.sum(PU_VComponents,axis = 0)))

    def Discounting(self,X):
        """Terms that pre-multiply the continution value in the value functions and envelope conditions.
        For the value functions these are just beta, but for the envelope conditions they are more complicated."""
        D = np.zeros((len(iVRange),X.shape[1] if X.ndim > 1 else 1))
        D[:] = beta

        return D

    def ForwardLooking(self,X,Xprime):
        '''iM, iVn, iV_a, iV_b, iV_c, iV, iYhat'''
        M  = (X[ih] * (X[iA]/mu - X[iw])/psi_1 + (1-upsilon)*Xprime[iM]**psi_2)**(1./psi_2)

        Vn = ( -self.log(self.b)-X[ih]**(1+gamma)/(1+gamma)+psy
                + beta * (1-upsilon)*(1-kappa/(1+kappa)*Xprime[iq]*Xprime[iM])*Xprime[iVn])

        if Vn.__class__ == np.ndarray and Vn[0].__class__ == np.float64  and any(np.isnan(Vn)):
            print any(np.isnan(X[ih]))
            print any(np.isnan(Xprime[iq]))
            print any(np.isnan(Xprime[iM]))
            print any(np.isnan(Xprime[iVn]))
            raise RuntimeError("Vn NaN")

        return np.vstack((M,
                    Vn,
                    self.PeriodUtility(X) + beta * Xprime[[iV_a,iV_b,iV_c,iV]],
                    X[iY]
                    ))


    def SolvePolicyIt(self,grid,Findices):
        '''Solve for value function with policy iteration idea'''
        X = np.zeros((self.nXYZ,grid.shape[1]))
        X[:self.nX] = grid
        X[self.nX:self.nXY] = self.F(X[self.interpstates])
        X[self.nXY:] = self.Static(X)

        BP = PI.PolyIntBasis(grid[self.interpstates],self.F.ord).T


        pdutil = self.PeriodUtility(X)
        # print pdutil

        EBPPrime = 0.0

        for quad_i in range(np.prod(self.quad['n'])):
            Stateprime = self.StateTrans(X,self.quad['x'][quad_i])
            EBPPrime += self.quad['w'][quad_i] * PI.PolyIntBasis(Stateprime[self.interpstates],self.F.ord).T


        coefs = np.zeros((pdutil.shape[0],BP.shape[1]))
        for i in range(pdutil.shape[0]):
            coefs[i] = np.linalg.lstsq(BP - beta * EBPPrime,pdutil[i],rcond=None)[0]


        self.F.coef[Findices] = coefs


    def Observation(self,X):
        """ Pick out the unemployment rate, GonY, log(Y), pi, h, M"""

        GDP = 1/X[iMU] + X[iG]
        EUrate = upsilon*(1-X[iq]*X[iM])

        return np.vstack(( X[iu],
                            X[iG]/GDP,
                            self.log(GDP),
                            X[ipi],
                            X[ih],
                            X[iM],
                            EUrate
                            ))

#################
# Numerical Options
#################
poly_ord = 2
F = PI.Poly(poly_ord)

# get quadrature grid
nquad = np.array((3,3))


T = 2500  # number of periods to simualte for
nGrid = 100  # number of points in the constructed grid

np.random.seed(8290372)

#-- initialization -- create a model object
VarNums = [5, 7, 8]
interpstates = WelfareStates
M = OptStabDynamics(VarNums,shock_covar,F,IncProc,b,tau,nquad,interpstates)
M.TestSteadyState()


def testfcnCycles(b,tau,M):

    M.Set_b_tau(b,tau)

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
    X = grid

    # add some Sobol points
    GridWidth = 1.5 * ( X.max(axis = 1) - X.min(axis = 1) )
    Sob = Drivers.SobolGrid(X.mean(axis = 1),GridWidth, nSobol = 20)
    X = np.hstack((X,Sob))


    # -- Solve for Welfare on the grid --
    #solve for welfare
    WelfareFindices = np.array((iV_a,iV_b,iV_c,iV))-M.nX
    M.SolvePolicyIt(grid,WelfareFindices)
    Welfare_Cyc = M.F(M.SteadyState()[M.interpstates][np.newaxis].T)[WelfareFindices].flatten()

    return Welfare_Cyc[3]

def testfcnSS(b,tau,M):
    M.Set_b_tau(b,tau)
    Welfare_SS = M.SteadyState()[iVRange].flatten()
    return Welfare_SS[3]


#-------------------------------------------------------------------
#---------------- Main script to do the optimization----------------
#--------------------------------------------------------------
if __name__ == '__main__':
    print "Steady State"
    res1 = optimize.minimize(lambda x: -testfcnSS(x[0],x[1],M), [0.75,0.35],method = 'Nelder-Mead' )
    print res1


    print "With Cycles"
    res2 = optimize.minimize(lambda x: -testfcnCycles(x[0],x[1],M), res1["x"],method = 'Nelder-Mead' )
    print res2


    np.savez(os.path.join(OutDir,'SearchResults_v' + str(VERSION) + '.npz'),OptSS = res1['x'], OptCyc = res2['x'])
