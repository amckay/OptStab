"""Model code for RBC model.

State transition:
K' = (1-delta) K  + Y - MU^(-1/gamma)
z' = rho z + eps

Static:
Y = exp(z) * K^alpha * L^(1-alpha)
R = exp(z) * alpha * K^(alpha-1) * L^(1-alpha) + 1 -delta
W = exp(z) * (1-alpha) * K^(alpha) * L^(-alpha)
L =   (W *MU)^(1/psi)

Forward looking:
MU = beta E_z' [ R' MU']
"""

import numpy as np
import ModelFramework.PolyInt as PI
import matplotlib.pyplot as plt
import ModelFramework.Drivers as Drivers
from ModelFramework.AutoDiff import exp, log
from ModelFramework.DetSimul import DetSimul

alpha = 0.3
beta = 0.99
delta = 0.02
gamma = 1.0
rho = 0.9
psi = 1.0


shock_covar = np.array([0.01**2])

#indexing
iK, iz, iMU, iY, iR, iW, iL = range(7)
VarNums = [2,1,4]

class RBC(Drivers.Model):




	def SteadyState(self):
		zstar = 0.0
		Rstar = 1./beta
		KL = ((Rstar - 1. + delta)/ ( alpha * exp(zstar)))**(1.0/(alpha - 1.0))
		Wstar = (1.0 - alpha) * exp(zstar) * KL**alpha

		YL = exp(zstar) * KL**alpha

		Lstar = ((YL - delta * KL) * Wstar**(-1./gamma))**(-1.0/(psi/gamma +1))
		Ystar = YL * Lstar
		Kstar = KL * Lstar
		Cstar = Ystar - delta * Kstar
		MUstar = Cstar**(-gamma)

		Xstar = np.array((Kstar,zstar,MUstar,Ystar,Rstar,Wstar,Lstar))
		

		return Xstar

	def StateTrans(self,X,shockprime):
		"""order of states is K, z"""
		return np.vstack( ( (1-delta)*X[iK] + X[iY] - X[iMU]**(-1.0/gamma),
							rho * X[iz] + shockprime[0]))

	def Static(self,X):
		"""order of static variables is Y, R, W, L"""

		W = (exp(X[iz]) * (1-alpha) * X[iK]**alpha * X[iMU]**(-alpha/psi))**(1./(1+alpha/psi))
		L =   (W *X[iMU])**(1.0/psi)
		Output = exp(X[iz]) * X[iK]**alpha * L**(1-alpha)
		R = alpha * Output/X[iK] + 1 - delta

		return  np.vstack(( Output, R, W, L))



	def ForwardLooking(self,X,Xprime):
		"""Returns marginal utility implied by the right-hand side of the Euler equation assuming that 
		the economy is in steady state and we know the forward looking variables in the next period."""
		
		
		MUprime = Xprime[iMU]
		Rprime = Xprime[iR]
		EMUprime = Rprime * MUprime 

		return np.array(beta * EMUprime)



# create a model object
nquad = 3
poly_ord = 2
M = RBC(VarNums,shock_covar,PI.Poly(poly_ord),nquad)
M.TestSteadyState()
PertSol = M.Perturb(SetF = True)


# Now we solve the model using the EDS grid
T = 1000
nGrid = 40
grid, Xsim = M.SolveEDS(T,nGrid)



plt.subplot(2,1,1)
plt.plot(Xsim[iMU]**(-1.0/gamma))
plt.title('Consumption')
plt.subplot(2,1,2)
plt.plot(Xsim[iL])
plt.title('Labor supply')
#plt.show()



DynareTestCase = [
	   np.array(((0.955551080144196,   2.932803087377220),
   		(0.0,   0.900000000000000))),
	   np.array(((0.,),(1.,))),
	np.array(((  -0.010518730616035,  -0.136867171144491),
	  ( 0.017857583827983 ,  3.483285050608325),
	  (-0.000991190985615 ,  0.041860302202998),
	  ( 0.026463258712741 ,  1.561228940597153),
	  (-0.006529187442374 ,  0.521823942066058)))   ]

for i in range(3):
	assert np.allclose(PertSol[i].flatten(),DynareTestCase[i].flatten()),  'failed for i = {}'.format(i)



print 'Deterministic Simulation'
T = 250
eps = np.zeros((1,T))
eps[0,1] = 0.002
X = M.SteadyState()
X0 =  np.tile(X,(T,1)).T

Xsim = M.DetSimul(X0,eps)


