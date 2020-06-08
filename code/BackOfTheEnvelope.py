import os
from OptStab import *
from scipy.optimize import bisect

OutDir = os.path.join(os.path.pardir,'data','results')

with np.load(os.path.join(OutDir,'SearchResults_v1.npz')) as X:
    bstar, taustar = X['OptCyc']

M.Set_b_tau(bstar,taustar)
Ystar = M.Ybar

def test_b(db):
    M.Set_b_tau(bstar+db,taustar)
    return M.Ybar - (Ystar/1.005)


def test_tau(dtau):
    M.Set_b_tau(bstar,taustar+dtau)
    return M.Ybar - (Ystar/1.005)


db = bisect(test_b,0,0.12)
b2 = bstar + db
print('found b2')

dtau = bisect(test_tau,0,0.05)
tau2 = taustar + dtau
print('found tau2')



def GetSTD():
    M.exp = adexp
    M.log = adlog
    PertSol = M.Perturb(SetF = True)
    M.exp = exp
    M.log = log

    SS = {'state0': M.SteadyState()[:M.nX][np.newaxis].T,'nrep':10}
    # grid, Xsim = M.SolveEDS(T/10,nGrid,damp = [0.995,0.995],verbose = True, PlotGrids = False, SimStart = SS)
    T = 500
    Xsim = M.Simulate(T)
    print('Done simulating')


    Q_u = np.zeros(T-1)
    Q_eps = np.zeros(T-1)
    for t in range(T-1):
        EUprime = upsilon*(1-Xsim[iq,t+1]*Xsim[iM,t+1])
        Q_u[t] = (1 - EUprime*(1-1.0/M.b(Xsim[iu,t+1])))
        Q_eps[t] = IncProc.get_Moments(Xsim[iu,t],M.ubar,M.tau,M.b(Xsim[iu,t]))[2]

    return np.std(Q_u), np.std(Q_eps), np.std(Q_u*Q_eps)




print("Star")
M.Set_b_tau(bstar,taustar)
print(GetSTD())

print("b2")
M.Set_b_tau(b2,taustar)
print(GetSTD())


print("tau2")
M.Set_b_tau(bstar,tau2)
print(GetSTD())
