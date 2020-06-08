import numpy as np
from OptStab_IncomeProcess import IncomeProcess
import json
import matplotlib.pyplot as plt

with open('../../data/calibration/Moments.txt','r') as infile:
    TargetMoments = json.load(infile)

delta = 1.0/200.0
ubar = TargetMoments['meanu']
ucalib = ubar

b = 0.81
tau = 0.151

IncProcVersion = 1
if IncProcVersion == 1:
    # This version uses x = u - ucalib
    IP_u2x = np.array((0.0, 16.73, 0.0))
else:
    # This version uses x = (1-u)/(1-ubar) - 1
    IP_u2x = np.array((0.0, -15.7596708786, 0.0))

IPP = {'sigma': (0.0403,0.0966,0.0966), 'p': (0.,0.00727,0.00727),
       'mu_0': (0.,0.266,-0.184), 'mu_x': (0.,1.,1.), 'neps': (5,7,7)}

IncProc = IncomeProcess(IPP,IP_u2x,ucalib,IncProcVersion)


#--
N = 100000
T = 2000

Sigmas = IncProc.sigma.flatten()
Means = IncProc.__get_MU__(0.0)[np.cumsum(np.hstack(((0.,),IncProc.neps)))[:-1].astype(int)].flatten()
i = np.zeros(N,int)
dead = np.zeros(N,bool)

log_var = np.zeros(T)

CEscale = np.ones(N)
CEscale[np.random.rand(N) < ubar] = b
CEscalevar = np.var(np.log(CEscale))


lal = np.zeros(N)

for t in range(T):
    i[:] =  np.random.choice([0,1,2], N, p=IncProc.p.flatten())
    lal += Means[i] + Sigmas[i] * np.random.randn(N)
    dead[:] = np.random.rand(N) < delta
    lal[dead] = 0.0
    log_var[t] = CEscalevar + (1-tau)**2 * np.var(lal)
    
plt.plot(log_var)    