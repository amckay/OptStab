print 'This is the version for low upsilon.'
execfile('OS_Version_1.py')
#
# import os
# OutDir = os.path.join(os.path.pardir,os.path.pardir,'data','results')
# with np.load(os.path.join(OutDir,'SearchResults_v1.npz')) as X:
#     bstar, taustar = X['OptSS']
#
#
# M.Set_b_tau(bstar,taustar)
#
# utarget = M.ubar
utarget = 0.06412234744773654


upsilon = 0.033448289852720
# upsilon  = utarget * qM /(1-qM - utarget * (1-qM))
qM = 1-utarget/(utarget + upsilon*(1-utarget))


# load
import os
OutDir = os.path.join(os.path.pardir,os.path.pardir,'data','results')
with np.load(os.path.join(OutDir,'SearchResults_v1.npz')) as X:
    bstar, taustar = X['OptSS']


def Calibw(X,Mode = 'Residuals'):
    Y, h = X # unpack
    A = S = 1.0  # steady state values


    Vn = -log(bstar) / (1-beta*(1-upsilon)*(1-kappa/(1+kappa) *qM) )
    q = (qM*Vn)**(1./(kappa+1))
    M = qM/q

    w = A/mu - upsilon*psi_1*M**psi_2/h
    H = 1-utarget-(1-upsilon)*(1-utarget)
    J = psi_1 * M**psi_2 * H
    h2 = ((1-taustar) * w * S  / A  * Y/ (Y-J))**(1.0/(1+gamma))
    Y2 = A * h2 * (1-utarget)

    if Mode.lower() == 'residuals':
        return np.array((Y2 - Y, h2 - h))
    elif Mode.lower() == 'values':
        return  w, Y, J, h, M, Vn
    else:
        raise Exception('Mode not recognized.')

sol = optimize.root(Calibw,np.array((1.0, 1.0)),jac=False)
wcalib, Ycalib, Jcalib, hbar, Mcalib, Vncalib = Calibw(sol.x,Mode = 'values')
print('upsilon = {0}'.format(upsilon))
print('wcalib = {0}'.format(wcalib))
