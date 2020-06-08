import os
from OptStab import *

OutDir = os.path.join(os.path.pardir,'data','results')

with np.load(os.path.join(OutDir,'SearchResults_v1.npz')) as X:
    OptCyc = X['OptCyc']

M.Set_b_tau(*OptCyc)

M.exp = adexp
M.log = adlog
PertSol = M.Perturb(SetF = True)
M.exp = exp
M.log = log

SS = {'state0': M.SteadyState()[:M.nX][np.newaxis].T,'nrep':10}
grid, Xsim = M.SolveEDS(T/10,nGrid,damp = [0.995,0.995],verbose = True, PlotGrids = False, SimStart = SS)

Mom = M.GetMoments()
print(np.sqrt(Mom['cov'][2,2]))
