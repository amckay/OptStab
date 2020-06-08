import numpy as np
import ModelFramework.Drivers as Drivers
import ModelFramework.EDS
import ModelFramework.PolyInt as PI
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
from Arguments import CLArgs


VERSION = CLArgs.ver



# indexing
(iEAlph, iA, iS, iMP, ieG, iulag, iElogAlpha,
 iMU, ipA, ipB, iJhat, iVn, idVndb, idVndtau, iV_a, iV_b, iV_c, iV,
 ipstar, ipi, iI, iG, iY, iSA, ih, iM, iw, iu, iJ, iH, iq) = range(31)

iVRange = np.array([ iV_a, iV_b, iV_c, iV])

class OptStabDynamics(Drivers.Model):


    def __init__(self,VarNums,shock_covar,F,IncProc,b,tau,nquad,interpstates):
        self.log = log
        self.exp = exp
        self.IP = IncProc
        self.Set_b_tau(b,tau)
        self.VarNums = VarNums


        # call the parent init
        super(OptStabDynamics,self).__init__(VarNums,shock_covar,F,nquad,interpstates)




    def Set_b_tau(self,bSS,tau,bCyc=0.0):

        self.bSS = bSS
        self.tau = tau
        self.bCyc = bCyc
        self.wbar = np.maximum(0.99*wcalib,wcalib * (bSS/bcalib)**WageElasticity)

        # # set dubar/db and d ubar / d tau
        # u0 = self.SteadyState_Labor_Market(bSS,tau)[4]
        # ub = self.SteadyState_Labor_Market(bSS+0.01,tau)[4]
        # utau= self.SteadyState_Labor_Market(bSS,tau+0.01)[4]
        # self.dubardb = (ub-u0) / 0.01
        # self.dubardtau = (utau-u0) / 0.01


        SS = self.SteadyState()
        self.Ibar = SS[iI]
        self.hbar = SS[ih]
        self.ubar = SS[iu]
        self.Ybar = SS[iY]
        self.Jbar = SS[iJ]
        self.Gbar = SS[iG]
        self.XSS = SS



    def b(self,u=None):
        if u is None:
            return self.bSS
        else:
            return np.maximum(np.minimum(self.bSS + self.bCyc * (u - self.ubar),0.95),0.05)


    def SteadyState_Labor_Market(self,b,tau):

        # six equations in h, Y, H, q, J, u
        # guess Y, J, and Vn
        def g(X,Mode = 'Residuals'):
            Y, J, Vn = X  #unpack
            h = ( (1-tau) * self.wbar * Y/ (  (Y-J)  )) **(1./(1.+gamma))
            M = ((1/mu  - self.wbar)*h/psi_1/upsilon)**(1/psi_2)
            q = (M*Vn)**(1.0/kappa)
            u = upsilon*(1 - q * M) / (q*M + upsilon*(1 - q * M) )
            H = 1-u - (1-upsilon)*(1-u)
            # check Y, J, Vn
            Y2 =   h * (1-u)
            J2 = psi_1 * M**psi_2 * H
            Vn2 = (-log(b) -   h**(1+gamma)/(1+gamma) + psy )/(1-beta*(1-upsilon)*(1-kappa/(1+kappa)*q*M) )
            if Mode.lower() == 'residuals':
                return np.array((Y2 - Y, J2-J, Vn2-Vn))
            elif Mode.lower() == 'values':
                return Y, J,  h, q, u, H, M, Vn
            else:
                raise Exception('Unrecognized value for Mode')


        sol = optimize.root(g,(Ycalib,Jcalib,Vncalib),jac=False)
        Y, J, h, q, u, H, M, Vn = g(sol.x,Mode = 'values')
        return Y, J, h, q, u, H, M, Vn

    def SteadyState(self,b = None, tau = None):
        if b is None:
            b = self.bSS

        if tau is None:
            tau = self.tau

        Y, J, h, q, u, H, M, Vn = self.SteadyState_Labor_Market(b,tau)

        w = self.wbar
        A = 1.
        MP_shock = 1.
        G_shock = 1.
        etam = 1.0
        S = 1.0
        SA = S
        pi = 1.
        pstar = 1.

        IP_Moments = self.IP.get_Moments(u,u,tau,b)
        ElogAlpha = (1-delta)*IP_Moments[1]/delta
        EAlph = float(delta / (1.0 - (1.0 - delta)*IP_Moments[0]))


        I = 1.0/  beta
        if not VERSION == 12:
            EU = upsilon*(1-q*M)
            I = I/((1-EU) + EU/b)
        if not VERSION == 13:
            I = I / IP_Moments[2]

        C = (Y-J)/(1+SSGRatio)
        G = SSGRatio * C
        pB = Y/(1-(1-theta)/I)
        pA = pB

        MU = 1.0/C
        dVndb = ( -1./b  ) / (1-beta*(1-upsilon)*(1-q*M))

        dVndtau =  (h**(1+gamma) / (1-tau) / (1+gamma)  )/ (1-beta*(1-upsilon)*(1-q*M))


        X = np.array((EAlph, A, S,  MP_shock, G_shock, u, ElogAlpha,
                        MU, pA, pB, J, Vn, dVndb, dVndtau,0.,0.,0.,0.,
                        pstar, pi, I, G, Y, SA, h, M, w, u, J, H, q ))


        V = np.dot(np.linalg.inv(np.eye(len(iVRange)) - np.diag(self.Discounting(X).flatten()) ) ,
                       self.PeriodUtility(X,SS=True)    )

        X[iVRange] = V.flatten()

        return X



    def StateTrans(self,X,shockprime):
        """Order of states is:
           0    1   2       3          4  5     6
        (EAlph, A, S,  MP_shock, G_shock, ulag, ElogAlpha   )"""

        if VERSION == 16:
            IP_Moments = self.IP.get_Moments(0.0*X[iu]+1.0*self.ubar,self.ubar,self.tau)
        elif VERSION == 31:
            IP_Moments = self.IP.get_Moments(1.25*(X[iu]-self.ubar)+1.0*self.ubar,self.ubar,self.tau)
        else:
            IP_Moments = self.IP.get_Moments(X[iu],self.ubar,self.tau,self.b(X[iu]))


        ElogEpsilonPrime = IP_Moments[1]

        return np.vstack((
            (1-delta) * X[iEAlph] * IP_Moments[0][0] + delta,
            self.exp(rhoA * self.log(X[iA]) + shockprime[0] ),
            X[iSA],
            self.exp(rhoMP * self.log(X[iMP]) + shockprime[1]),
            self.exp(rhoG * self.log(X[ieG]) + shockprime[2]),
            X[iu],
            (1-delta) * (X[iElogAlpha] + ElogEpsilonPrime)
            ))


    def Static(self,X):
        """order of static variables is
          0     1   2   3  4   5  6  7  8  9  10 11  12
        (pstar, pi, I, G, Y, SA, h, M, w, u, Mf, H, q)"""


        pstar = np.maximum(X[ipA]/X[ipB],0.0001)
        pi = np.maximum((1-theta)/(1-theta * pstar**(1/(1-mu))),0.001)**(1-mu)

        C = np.maximum(1/X[iMU],0.001)

        if VERSION == 10 or VERSION == 32:
            G = np.maximum(self.Gbar * np.ones(X.shape[1])* X[ieG],0.001)
        else:
            G = np.maximum(SSGRatio * C * X[ieG],0.001)

        if VERSION == 4 or  VERSION == 15:
            SACap = 1.1
        else:
            SACap = 1.03

        SA = np.minimum((1-theta) * X[iS] * pi**(-mu/(1-mu)) + theta * pstar**(mu/(1-mu)),SACap)
        if SA.__class__ == np.ndarray and (SA[0].__class__ == np.float64 or SA[0].__class__ == np.ndarray ) and any(np.isnan(SA)):
            print("Found nan in SA")
            print(X[iS])
            print(pi)
            print(pstar)
        Jtemp = np.maximum(X[iJhat],0.00005*C)
        Y = C + G + Jtemp
        # if  np.any((SA * Y * self.hbar**Phi) /  (1-self.ubar)  < 0.0):
        #     print 'weirdness'
        #     print Y

        # Here we are imposing a wage rule although it is not obvious because
        # we have done some substitution.  We use the wage rule and the labor supply
        # rule to solve for h in terms of 1-u.  We then substitute in for h in
        # the production function so we have Y as a function of u (and some other things
        # but not h) and then we solve for u in terms of Y.
        #
        # Derivation uses equations
        # h^\gamma &= \frac{(1-\tau) w}{A h (1-J/Y)} \\
        # w &= \bar w \eta_t^A \left( \frac{1-u}{1-\bar u} \right)^\zeta \\
        # Y &= A h (1-u) \\
        # A &= \eta^A / S \\
        # steps
        # h &= \left(\frac{(1-\tau) w}{A  (1-J/Y)} \right)^{1/(1+\gamma)} \\
        # Y &= \frac{\eta }{ S} \left(\frac{(1-\tau) w}{A  (1-J/Y)} \right)^{1/(1+\gamma)} (1-u) \\
        # Y &= \frac{\eta }{ S} \left(\frac{(1-\tau) \bar w \eta^{\zeta_1} }{A  (1-J/Y) \left(1-\bar u\right)^{\zeta_2} }\right)^{1/(1+\gamma)} (1-u)^{1+\frac{\zeta_2}{1+\gamma}} \\
        # 1-u &= \left[ \frac{Y}{A} \left(\frac{(1-\tau) \bar w \eta^{\zeta_1} }{A  (1-J/Y) \left(1-\bar u\right)^{\zeta_2} }\right)^{-1/(1+\gamma)} \right]^{\frac{1+\gamma}{1+\gamma+\zeta_2}}

        u = 1- (Y/(X[iA]/SA)*(((1-self.tau)*self.wbar*X[iA]**zeta_1)/((X[iA]/SA) *(1-Jtemp/Y) * (1-self.ubar)**zeta_2))**(-1./(1.+gamma)) )**((1.+gamma)/(1.+gamma+zeta_2))

        w = self.wbar * X[iA]**zeta_1 * ((1-u)/(1-self.ubar))**zeta_2

        h = ((1-self.tau)*w/((X[iA]/SA) * (1-Jtemp/Y) ))**(1./(1+gamma))

        if h.__class__ == np.ndarray and (h[0].__class__ == np.float64 or h[0].__class__ == np.ndarray ) and any(np.isnan(h)):
            print("Found nan in h")
            print(Y)
            print(C)
            print(G)
            print(X[iJhat])
            print(X[iA])
            print(SA)
            raise RuntimeError("Found nan in h")
        # We use the subsitutions above to recover u from h
        numsearchers = X[iulag]+upsilon*(1-X[iulag])
        u = np.minimum(u,numsearchers-0.005)

        # Use:
        # q^kappa = M V^n    and   u  = (ulag + upsilon (1-ulag))*(1-qM)
        # to solve for M given u and Vn
        M = ( np.maximum(numsearchers - u,0.001)/ np.maximum(numsearchers,0.001)  * np.maximum(X[iVn] , 0.0001 )**(-1./kappa)  )**(kappa/(1.0 + kappa))
        # use   u  = (ulag + upsilon (1-ulag))*(1-qM)  to find q
        q = (1- u / numsearchers)/M
        if any(q < 0):
            print("----")
            print("Y")
            print(Y)
            print("A")
            print(X[iA])
            print("SA")
            print(SA)
            print("h")
            print(h)
            print("q")
            print(q)
            print("u")
            print(u)
            print("numsearchers")
            print(numsearchers)
            print("M")
            print(M)
            raise RuntimeError("negative q")
        H = numsearchers  - u
        I = self.Ibar * pi**omegapi * X[iMP] * ((1-u)/(1-self.ubar))**(omegau)
        J = psi_1 * M**psi_2 * H

        return np.vstack((pstar, pi, I, G, Y, SA, h, M, w, u, J, H, q))



    def ForwardLooking(self,X,Xprime):
        """Returns marginal utility implied by the right-hand side of the Euler equation assuming that
        the economy is in steady state and we know the forward looking variables in the next period."""

        EUprime = upsilon*(1-Xprime[iq]*Xprime[iM])

        if VERSION == 12:
            Q_u = 1.0
        else:
            Q_u = (1 - EUprime*(1-1.0/self.b(Xprime[iu])))

        if VERSION == 13:
            Q_eps = 1.0
        elif VERSION == 16:
            Q_eps = IncProc.get_Moments(0.0*X[iu]+1.0*self.ubar,self.ubar,self.tau)[2]
        elif VERSION == 31:
            Q_eps = IncProc.get_Moments(1.5*(X[iu]-self.ubar)+1.0*self.ubar,self.ubar,self.tau)[2]
        else:
            Q_eps = IncProc.get_Moments(X[iu],self.ubar,self.tau,self.b(X[iu]))[2]




        MU =  (     beta *  X[iI]/Xprime[ipi] * Xprime[iMU]  * (1 - Xprime[iu]*(1-self.b(Xprime[iu])))
                    / (1-X[iu]*(1-self.b(X[iu]))) * Xprime[iEAlph]/X[iEAlph] * Q_u * Q_eps )

        pA = Xprime[ipi]**(-mu/(1-mu)) * (Xprime[ipi]/Xprime[iI]) * Xprime[ipA]
        pB = Xprime[ipi]**(-mu/(1-mu)) * (Xprime[ipi]/Xprime[iI]) * Xprime[ipB]



        pA = mu * X[iY] * (X[iw]/X[iA] + psi_1 * X[iM]**psi_2/X[ih]/X[iA] - (1-upsilon)*psi_1*Xprime[iM]**psi_2/X[ih]/X[iA] )  + (1-theta) * pA
        pB = X[iY] + (1-theta) * pB
        Jhat = X[iJ]
        Vn = ( -self.log(self.b(X[iu]))-X[ih]**(1+gamma)/(1+gamma)+psy
                + beta * (1-upsilon)*(1-kappa/(1+kappa)*Xprime[iq]*Xprime[iM])*Xprime[iVn])

        dudq = -X[iM]*(X[iulag]+upsilon*(1-X[iulag]))
        dhdq = -zeta_2 / (1+gamma) * X[ih]/(1-X[iu]) * dudq
        dhdubar =  zeta_2 / (1+gamma) * X[ih]/(1-self.ubar)

        dubardq = -self.XSS[iM]*(self.XSS[iulag]+upsilon*(1-self.XSS[iulag]))

        dqbardb = (1./kappa)*self.XSS[iq]/self.XSS[iVn]*self.XSS[idVndb]
        dubardb = dubardq * dqbardb

        dqbardtau = (1./kappa)*self.XSS[iq]/self.XSS[iVn]*self.XSS[idVndtau]
        dubardtau = dubardq * dqbardtau

        dVndb = ( (1+X[ih]**gamma * dhdq / kappa * X[iq] / Vn)**(-1)
                    *(-1./self.b(X[iu]) - X[ih]**gamma*dhdubar* dubardb  + beta * (1-upsilon) *(1-Xprime[iq]*Xprime[iM])*Xprime[idVndb])
                )
        dVndtau = ( (1+X[ih]**gamma / kappa * X[iq] / Vn * dhdq)**(-1)*
                    (   X[ih]**(1+gamma) / (1-self.tau) / (1+gamma)
                        -X[ih]**gamma * dhdubar * dubardtau
                        + beta * (1-upsilon) *(1-Xprime[iq]*Xprime[iM])*Xprime[idVndtau]
                    )
                   )

        return np.vstack((MU,pA,pB,Jhat,Vn,dVndb,dVndtau,self.PeriodUtility(X) + self.Discounting(X) * Xprime[iVRange] ))

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

    def PeriodUtility(self,X, SS = False):
        """ Period contributions to welfare and envelope conditions for
        iVn, iV_a, iV_b, iV_c, iV, iEnv_EAlph, iEnv_ElogAlpha, iEnv_S, iEnv_ulag_Hos, iEnv_ulag_Risk"""


        if SS:
            b = self.bSS
        else:
            b = self.b(X[iu])

        if np.any((-self.log(b)- gamma_0 * (X[ih]**(1+gamma)/(1+gamma) -psy)).ravel() ) < 0:
            warnings.warn('Negative search effort encountered. A')

        if np.any(X[iq].ravel() ) < 0:
            warnings.warn('Negative search effort encountered. B')

        if X[iq].__class__ == np.ndarray and (X[iq][0].__class__ == np.float64 or X[iq][0].__class__ == np.ndarray ) and any(np.isnan(X[iq])):
            warnings.warn('NaN search effort encountered.')

        if np.any(X[ih].ravel() ) < 0:
            warnings.warn('Negative work effort encountered.')


        numsearchers = upsilon*(1-X[iulag]) + X[iulag]

        PU_VComponents = np.vstack( (  (1-self.tau)*X[iElogAlpha] - self.log(X[iEAlph]),
                                    X[iu] * self.log(b) - self.log(1- X[iu]*(1-b)),
                                    -self.log(np.maximum(X[iMU],0.001)) - (1-X[iu])*gamma_0 * (X[ih]**(1+gamma)/(1+gamma)-psy) - numsearchers * X[iq]**(1+kappa)/(1+kappa) + chi * self.log(np.maximum(X[iG],0.001))
                                    ))


        return np.vstack((PU_VComponents, np.sum(PU_VComponents,axis = 0)))

    def Discounting(self,X):
        """Terms that pre-multiply the continution value in the value functions and envelope conditions.
        For the value functions these are just beta, but for the envelope conditions they are more complicated."""
        D = np.zeros((len(iVRange),X.shape[1] if X.ndim > 1 else 1))
        D[:] = beta

        return D


    def SolvePolicyIt(self,grid,Findices):
        '''Solve for value function with policy iteration idea

        Note that this works for the value function but not the envelope conidtions
        because the discounting is just beta.
        '''
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




execfile('OptStab_calibration.py')


#################
# Numerical Options
#################
poly_ord = 2
F = PI.Poly(poly_ord)

# get quadrature grid
nquad = np.array((3,3,3))


T = 2500  # number of periods to simualte for
nGrid = 100  # number of points in the constructed grid

np.random.seed(8290372)


#-- initialization --

# create a model object
VarNums = [7, 11, 13]
interpstates = DynamicsStates
M = OptStabDynamics(VarNums,shock_covar,F,IncProc,b,tau,nquad,interpstates)
M.TestSteadyState()
