"""Functions based on the CompEcon Matlab toolbox written by Miranda and Fackler."""



import numpy as np

def qnwnorm(n,mu,var):
    """QNWNORM Computes nodes and weights for multivariate normal distribution

    [x,w] = qnwnorm(n,mu,var);
    INPUTS
    n   : len d numpy array of number of nodes for each variable
    mu  : len d numpy array with means
    var : d by d numpy array with positive definite covariance matrix
    OUTPUTS
    x   : prod(n) by d numpy array of evaluation nodes
    w   : prod(n) by 1 numpy array of probabilities

    To compute expectation of f(x), where x is N(mu,var), write a
    function f that returns m-vector of values when passed an m by d
    matrix, and write [x,w]=qnwnorm(n,mu,var); E[f]=w'*f(x);

    USES: gridmake

    Adapted by Alisdair McKay from code of Miranda and Fackler
    Saturday, July 26, 2014 """

    if type(n) is int:
        d = 1
    else:
        d = n.size

    x = []
    w = []

    if d > 1:
        for i in range(d):
            xi,wi = _qnwnorm1(n[i])
            x.append(xi)
            w.append(wi)                
        
        wk = w[d-1]
    
        for i in range(d-2,-1,-1): 
            wk = np.kron(wk,w[i])
        w = wk

        x = _gridmake(x)

        x = np.dot(x,np.linalg.cholesky(var).transpose()) + np.kron(np.ones((np.prod(n),1)),mu)
                
    else:
        x,w = _qnwnorm1(n)
        x = np.sqrt(var) * x + mu

        x = x.flatten()
        w = w.flatten()        
  


    return (x,w)


def _gridmake(x):
    #recursively apply gridmake2 to all elements of list x
    nx = len(x)
    if nx == 1:
        return x
    elif nx == 2:
        return _gridmake2(x[0],x[1])
    else:
        return _gridmake([_gridmake2(x[0],x[1])] + x[2:])
    
def _gridmake2(x,y):
    """gridmake2 returns all combinations of rows of x and y into a single array with mx*my rows and nx + ny columns"""
    mx,nx = x.shape
    my,ny = y.shape
    
    return np.concatenate((np.kron(np.ones((my,1)),x), np.kron(y,np.ones((mx,1)))),axis = 1)




def _qnwnorm1(n):
    """QNWNORM1 Computes nodes and weights for the univariate standard normal distribution
    USAGE
    [x,w] = qnwnorm1(n);
    INPUTS
    n   : number of nodes
    OUTPUTS
    x   : n by 1 vector of evaluation nodes
    w   : n by 1 vector of probabilities
 

    Based on an algorithm in W.H. Press, S.A. Teukolsky, W.T. Vetterling
    and B.P. Flannery, "Numerical Recipes in FORTRAN", 2nd ed.  Cambridge
    University Press, 1992."""


    
    maxit = 100
    pim4 = 1/np.pi**0.25
    m = int(np.floor((n+1)/2))
    x = np.zeros((n,1))
    w = np.zeros((n,1))
    for i in range(m):
        # Reasonable starting values
        if i==0:
            z = np.sqrt(2*n+1)-1.85575*((2*n+1)**(-1./6.))
        elif i ==1:
            z = z-1.14*(n**0.426)/z
        elif i ==2:
              z = 1.86*z+0.86*x[0]
        elif i ==3:
              z = 1.91*z+0.91*x[1]
        else:
              z = 2*z+x[i-2]

    # root finding iterations
        its=0
        while its<maxit:
            its = its+1
            p1 = pim4
            p2 = 0
            for j in range(n):
                p3 = p2
                p2 = p1
                p1 = z*np.sqrt(2.0/(j+1.0))*p2-np.sqrt(j/(j+1.0))*p3

            pp = np.sqrt(2*n)*p2
            z1 = z
            z  = z1-p1/pp
            if abs(z-z1)<1e-14:
                break

        if its>=maxit:
            raise RuntimeError('failure to converge in _qnwnorm1')

        x[n-1-i] = z
        x[i] = -z
        w[i] = 2/(pp*pp)
        w[n-1-i] = w[i]

    w = w/np.sqrt(np.pi)
    x = x*np.sqrt(2)

    return (x,w)

