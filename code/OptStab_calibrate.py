from OptStab import *
from scipy.linalg import solve_discrete_lyapunov as lyapunov
from ModelFramework.AutoDiff import DoJac

DoSearch = False



# ------------------------------------------------------------------------
# This initial block calibrates the aggregate shocks (there are three)
# based on a first-order perturbation solution and three moments:
# 1) Standard deviation of u
# 2) Standard deviation of G/Y
# 3) Technology and monetary shocks contribute equally to var(u)
# ------------------------------------------------------------------------


M.exp = adexp
M.log = adlog
PertSol = M.Perturb(SetF = True)
LinearObs = DoJac(M.Observation,M.SteadyState())[0]
M.exp = exp
M.log = log



def GetVariances(shock_var):
    '''Returns the covariance matrix for the variables in the observation equation.
    Inputs: an array of variances for the aggregate shocks.
    Outputs: a 2D array covariance matrix corresponding to M.Observation
    Method: uses linearized solution.
    '''
    A, B,C =  PertSol
    C = np.vstack((np.eye(M.VarNums[0]), C  ))
    V0 = lyapunov(A, np.dot(np.dot(B , np.diag(shock_var)) ,B.T))
    V0 = np.dot(np.dot(C,V0),C.T)
    V0 = np.dot(np.dot(LinearObs,V0),LinearObs.T)
    return V0


# set up a system of equations based on moments for aggregate variances
# and variance decomposition
A = np.zeros((3,3))
for i in range(3):
    tmp = np.zeros(3)
    tmp[i] = 1.0
    V0 = GetVariances(tmp)
    # var(u), var(G/Y)
    A[:,i] = np.diag(V0[:3,:3])


B = np.zeros((3,1))
B[0] = TargetMoments['StDevu']**2
B[1] = TargetMoments['StDevGY']**2

# here we impose that the technology shock is calibrated externally
A[2,:] = np.array((1.,0.,0.0))
B[2] = 0.0046**2


M.shock_covar = np.diag(np.linalg.solve(A,B).flatten())
ScaleFactor = np.sqrt(M.shock_covar[1,1] / M.shock_covar[0,0])



# ------------------------------------------------------------------------
# Next we have code to compare the model non-linear to a range of targets
# it is possible to search for parameters to fit the non-linear model
# to the targets, but for now I am not doing this.
# ------------------------------------------------------------------------

# calibrate the shock covariance matrix
def CheckCalibration(M,DoPlots = False):


    SS = {'state0': M.SteadyState()[:M.nX][np.newaxis].T,'nrep':10}
    grid, Xsim = M.SolveEDS(T/10,nGrid,damp = [0.995,0.995] , PlotGrids = False, SimStart = SS)


    Mom = M.GetMoments()
    print 'means'
    print Mom['mn']


    print 'st deviations'
    print '  Model       Target'
    ModSTD =  np.sqrt(np.diag(Mom['cov']))
    TarSTD = [TargetMoments['StDev'+series] for series in ['u', 'GY', 'Y', 'infl', 'HoursPerWorkerIndex']]

    names = ['Unemp.   ', 'G / Y.   ', 'log Y    ', 'infla.   ', 'log hrs. ']

    for x in zip(names,ModSTD,TarSTD):
        print x[0] + "%8.4g    %8.4g" % x[1:]


    print 'correlation (Y,pi) = {}'.format(Mom['cov'][2,3]/np.sqrt(Mom['cov'][2,2] * Mom['cov'][3,3]))



    if DoPlots:
        Y = M.Observation(Xsim)
        for i in range(len(names)):
            plt.subplot(2,3,i+1)
            plt.plot(Y[i])
            plt.title(names[i])

        plt.show()


    return np.array((ModSTD[:2]-TarSTD[:2]))

def CheckCalibrationWrap(x):
    tmp = np.array([x[0],ScaleFactor*x[0],x[1]])
    M.shock_covar = np.diag(tmp)**2/1e4
    print M.shock_covar
    return CheckCalibration(M,DoPlots = False)




if DoSearch:
    print("Now searching")
    x0 = np.sqrt(np.diag(M.shock_covar)[[0,2]]*1e4)
    res = optimize.broyden1(CheckCalibrationWrap,x0,f_tol=1e-4,alpha = -1./0.04)
    tmp = np.array([res[0],ScaleFactor*res[0],res[1]])
    M.shock_covar = np.diag(tmp)**2/1e4

CheckCalibration(M,DoPlots = False)

print 'Shock variances (x 1000):'
print M.shock_covar*1000
