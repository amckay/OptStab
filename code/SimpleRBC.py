iK, iz, iM, iY, iR = range(5)
VarNums = [2,1,2]

import numpy as np
delta, rho, beta, alpha, gamma  = 0.02, 0.9, 0.99, 0.3, 2.0
shock_covar = np.array([[0.001]])  # note it is a 2-D array


from ModelFramework.Drivers import Model

class RBC(Model):
	def StateTrans(self,X,epsprime):
		Kprime = (1 - delta) * X[iK] + X[iY] - X[iM]**(-1/gamma)
		zprime = rho * X[iz] + epsprime[0]
		return np.vstack((Kprime,zprime))

	def ForwardLooking(self,X,Xprime):
		return np.vstack((beta * Xprime[iR] * Xprime[iM],))	# Note the trailing comma.  (1.0) is a float, (1.0,) is a tuple.
															# This is necessary becase we only have one return variable
															# and we want to force the returned object to be a 2D array.

	def Static(self,X):
		Y = exp(X[iz]) * X[iK]**alpha
		R = alpha * exp(X[iz]) * Y / X[iK] + 1 - delta
		return np.vstack((Y,R))

	def SteadyState(self):
		z = 0
		R = 1/beta
		K = ( (R - 1 + delta)/alpha)**(1.0/(alpha-1))
		Y = K**alpha
		M = ( Y - delta * K)**(-gamma)
		return np.array((K,z,M,Y,R)) # note the order matches our indices iK, iz, etc.

	def Observation(self,X):
		'''Pick out Y and M and transform M to C'''
		return np.vstack(( X[iY], X[iM]**(-1/gamma)) )


nquad = 5

import ModelFramework.PolyInt as PI
mypoly = PI.Poly(3)

M = RBC(VarNums,shock_covar,mypoly,nquad)


from ModelFramework.AutoDiff import exp, log
M.TestSteadyState()
PertSol = M.Perturb(SetF = True)

T = 1000  # number of periods to simulate for
nGrid = 40 # number of grid points to use in the projection method
grid, Xsim = M.SolveEDS(T,nGrid)

Moments = M.GetMoments(T)
print 'Means'
print 'Y             C'
print Moments['mn']