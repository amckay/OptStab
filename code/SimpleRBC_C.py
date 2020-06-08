iK, iz, iC, iY, iR = range(5)
VarNums = [2,1,2]

import numpy as np
delta, rho, beta, alpha, gamma  = 0.02, 0.9, 0.99, 0.3, 2.0
shock_covar = np.array([[0.001]])  # note it is a 2-D array


from ModelFramework.Drivers import Model


class RBC(Model):
    def StateTrans(self,X,epsprime):
        Kprime = (1 - delta) * X[iK] + X[iY] - X[iC]
        zprime = rho * X[iz] + epsprime[0]
        return np.vstack((Kprime,zprime))

    def Static(self,X):
        Y = exp(X[iz]) * X[iK]**alpha
        R = alpha * exp(X[iz]) * Y / X[iK] + 1.0 - delta
        return np.vstack((Y,R))

    def SteadyState(self):
        z = 0.0
        R = 1/beta
        K = ( (R - 1 + delta)/alpha)**(1.0/(alpha-1))
        Y = K**alpha
        C = Y - delta * K
        return np.array((K,z,C,Y,R)) # note the order matches our indices iK, iz, etc.

    def ForwardLooking(self,X,Xprime):
        return np.vstack(((beta * Xprime[iR])**(-1.0/gamma) * Xprime[iC],)) 


    def ForwardLookingExpectation(self,State):
        X = np.vstack((State, self.F(State)))
        X = np.vstack((X, self.Static(X)))

        MUimplied = 0.0

        for quad_i in range(np.prod(nquad)):
            Xprime = self.StateTrans(X,quad_x[quad_i])
            Xprime = np.vstack((Xprime, self.F(Xprime)))
            Xprime = np.vstack((Xprime,self.Static(Xprime)))
            
            MUimplied += quad_w[quad_i] * (beta * Xprime[iR] * Xprime[iC]**(-gamma))

        Cimplied = MUimplied**(-1.0/gamma)

        return  Cimplied

    

nquad = 5
from ModelFramework.CompEcon import qnwnorm
quad_x, quad_w = qnwnorm(nquad,0.0,shock_covar)
if quad_x.ndim == 1:
            quad_x = quad_x[np.newaxis].T

import ModelFramework.PolyInt as PI
mypoly = PI.Poly(3)
M = RBC(VarNums,shock_covar,mypoly,nquad)

from ModelFramework.AutoDiff import exp, log
M.TestSteadyState()
PertSol = M.Perturb(SetF = True)
grid, Xsim = M.SolveEDS(1000,40)