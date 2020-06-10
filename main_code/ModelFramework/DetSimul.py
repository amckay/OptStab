import numpy as np
from numpy import zeros, max, abs
import warnings
from AutoDiff import DoJac



def DetSimul(X0,StateTrans,ForwardLooking,Static, eps, nX,nY,nZ, maxit, Dampening = 1, DampeningThresh = 0):
    """ Deterministic solution to the model represented by StateTrans,ForwardLooking,Static.

    
    X0 -- initial guess of solution (nx x T)  
          Order of variables is state, forwardlooking, static
    eps -- exogenous variables (neps x T)
    StateTrans,ForwardLooking,Static -- function that give model equations
    maxit -- maximum number of iterations to try
    dampening -- initial scale factor of update: X' = X + dampening * dX
    DampeningThresh -- residual threshold at which dampening is turned off.

    X  = DetSimul( X0, eps, fcn)

    X is nx x T  solution with same ordering as X0

    This program implements the algorithm described in Juillard (1996) 
    "DYNARE: A program for the resolution and simulation of dynamic models 
    with forwardd variables through the use of a relaxation algorithm." 

    Alisdair McKay
    Saturday, August 1, 2015
    """

    nx0 = X0.shape[0]
    # print 'nx0 = {}'.format(nx0)

    def eqmwrap(X,eps):

        def eqm_f(XX):

            
            # Unpack
            Xlag = XX[:nx0][np.newaxis].T
            Xcur = XX[nx0:2*nx0][np.newaxis].T
            Xprm = XX[2*nx0:][np.newaxis].T

            Y = Xcur - np.vstack((
                            StateTrans(Xlag,eps),                
                            ForwardLooking(Xcur,Xprm),
                            Static(Xcur)
                        ))


            return Y

        
        J, Y = DoJac(eqm_f,np.hstack([X[:,i] for i in range(3)]))
        


        Jlag = J[:,:nx0]
        Jcur = J[:,nx0:2*nx0]
        Jprm = J[:,2*nx0:]


        return Y, Jlag, Jcur, Jprm

    return __Solve__( X0, eps, eqmwrap, maxit, Dampening, DampeningThresh )

def __Solve__( X0, eps, eqmwrap, maxit, Dampening, DampeningThresh):
    """ Deterministic solution to the model represented by eqmwrap.

    
    X0 -- initial guess of solution (nx x T)
    eps -- exogenous variables (neps x T)
    eqmwrap -- function that gives residuals of model equations
    maxit -- maximum number of iterations to try
    dampening -- initial scale factor of update.
    DampeningThresh -- residual threshold at which dampening is turned off.

    X  = DetSimul( X0, eps, fcn)

    X is nx x T

    This program implements the algorithm described in Juillard (1996) 
    "DYNARE: A program for the resolution and simulation of dynamic models 
    with forwardd variables through the use of a relaxation algorithm." 

    Alisdair McKay
    Saturday, August 1, 2015
    """



    #initializations
    X = X0

    nx, T = X.shape
    C = zeros((nx,nx,T))
    d = zeros((nx,T))
    dX = zeros((nx,T))
    Fx = zeros((nx,T))


    for it in range(maxit+1):

        C[:,:,0] = zeros(nx)
        d[:,0] = zeros(nx)
        for t in range(1,T-1):
            Fx[:,[t]],S_L, S, S_P = eqmwrap(X[:,t-1:t+2],eps[:,t])
            C[:,:,t] = np.linalg.solve( S - np.dot(S_L,C[:,:,t-1]),S_P)
            # print Fx[:,[t]].shape
            # print np.dot(S_L,d[:,[t-1]]).shape
            # print S.shape 
            # print np.dot(S_L,C[:,:,t-1]).shape
            # print (S - np.dot(S_L,C[:,:,t-1]) ).shape
            d[:,t] = -np.linalg.solve(S - np.dot(S_L,C[:,:,t-1]) ,Fx[:,t] + np.dot(S_L,d[:,t-1]))
            
        
        
        residual = max(abs(Fx.flatten()))
        if residual < 1e-8:
            print 'done'
            break
        else:
            print 'iteration ' + str(it) + '...residual ' + str(residual)
            if residual < DampeningThresh:
                Dampening = 1
            
            

        
        
        dX[:,T-1] = zeros(nx)
        for t in range(T-2,0,-1):
            dX[:,t] = d[:,t] - np.dot(C[:,:,t],dX[:,t+1])
        
        
        
        X = X + Dampening * dX
        

    if it >= maxit:
        warnings.warn('simulation did not converge after {} iterations.'.format(maxit))
    


    return X

