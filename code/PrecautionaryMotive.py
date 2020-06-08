from OptStab import M, adexp, adlog, iu, iY, IncProc
import numpy as np
from numpy import exp, log
import os

VERSION = 1
OutDir = os.path.join(os.path.pardir,os.path.pardir,'data','results')
OutputDeviation = -0.001

with np.load(os.path.join(OutDir,'SearchResults_v'+str(VERSION)+'.npz')) as X:
    OptCyc = X['OptCyc']
    OptSS = X['OptSS']

M.Set_b_tau(*OptSS)

M.exp = adexp
M.log = adlog
PertSol = M.Perturb(SetF = True)
M.exp = exp
M.log = log

grid, Xsim = M.SolveEDS(2000,100,tol = 1e-5, maxit = 5000, damp = [0.96, 0.98], PlotGrids = False,verbose = True)


# find the deviation in policy that leads to a 0.1% decline in output
PolicyDeviations = OptSS.copy()
from scipy.optimize import newton
Ybase = M.SteadyState()[iY]
def checkb(b):
    M.Set_b_tau(b,OptSS[1])
    return (M.SteadyState()[iY]/Ybase-1 -OutputDeviation)

def checktau(tau):
    M.Set_b_tau(OptSS[0],tau)
    return (M.SteadyState()[iY]/Ybase-1 -OutputDeviation)


PolicyDeviations[0] = newton(checkb,OptSS[0] )
PolicyDeviations[1] = newton(checktau,OptSS[1] )
#-----------------
# Compute the volatility of the precautionary motive with different policies
Q = np.zeros((2,Xsim.shape[1]-1))
Q2 = Q.copy()
for t in range(Xsim.shape[1]-1):
    Q[0,t] = (1 - Xsim[iu,t+1]*(1-1.0/OptSS[0]))
    Q[1,t] = IncProc.get_Moments(Xsim[iu,t],M.ubar,OptSS[1])[2]
    Q2[0,t] = (1 - Xsim[iu,t+1]*(1-1.0/PolicyDeviations[0]))
    Q2[1,t] = IncProc.get_Moments(Xsim[iu,t],M.ubar,PolicyDeviations[1])[2]


stdOfPrecMotive = Q2.std(axis = 1)/Q.std(axis = 1) -1

print('Change in policy leads to percent decline in precautionary motive:')
print('b: {0}'.format(stdOfPrecMotive[0]))
print('tau: {0}'.format(stdOfPrecMotive[1]))
