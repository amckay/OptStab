import numpy as np
from scipy.linalg import qz
from numpy import mat, c_, r_, where, sqrt, newaxis
from numpy.linalg import solve
from numpy.matlib import diag
from AutoDiff import DoJac

def Linearize(f,g,s,x,ne):
    """Linearizes the equations at the steady state.

    The original system is assumed to be in the the form:
    .. math::
        E_t f(s_t,x_t,s_{t+1},x_{t+1})
        s_t = g(s_{t-1},x_{t-1},  \\epsilon_t)

    Inputs:
    f and g are functions as defined above that take column vector inputs and return column vectors.
    s and x are column vectors of steady state values for s and x
    ne is number of epsilon shocks.

    Outputs: f1, g1  are matrices of first derivatives.
    f1 is shape nx by 2*(ns + nx)
    the first ns columns are derivatives w.r.t. :math:`s_t',
    the next nx columns are derivatives w.r.t. :math:`x_t',
    and so on.

    g1 is shape ns by (ns + nx + ne)
    the first ns columns are derivatives w.r.t. :math:`s_t',
    the next nx columns are derivatives w.r.t. :math:`x_t',
    the last ne columns are derivatives w.r.t. :math:`\\epsilon_t',

    """
    ns = len(s)
    nx = len(x)
    nv = ns + nx

    def fwrap(X):
        return f(X[:ns],X[ns:nv],X[nv:nv+ns],X[nv+ns:2*nv])

    def gwrap(X):
        return g(X[:ns],X[ns:nv],X[2*nv:])

    X = np.vstack((s,x,s,x,np.zeros((ne,1))))

    f1 = DoJac(fwrap,X)[0]
    g1 = DoJac(gwrap,X)[0]


    # print 'dY_3 / dX_1 = {}'.format(f1[2,0])

    f1 = f1[:,:2*nv]
    g1 = np.hstack((g1[:,:nv], g1[:,2*nv:]))

    return f1, g1





def qzswitch(i, A2, B2, Q, Z):
    #print i, A2, B2, Q, Z
    Aout = A2.copy(); Bout = B2.copy(); Qout = Q.copy(); Zout = Z.copy()
    ix = i-1    # from 1-based to 0-based indexing...
    # use all 1x1-matrices for convenient conjugate-transpose even if real:
    a = mat(A2[ix, ix]); d = mat(B2[ix, ix]); b = mat(A2[ix, ix+1]);
    e = mat(B2[ix, ix+1]); c = mat(A2[ix+1, ix+1]); f = mat(B2[ix+1, ix+1])
    wz = c_[c*e - f*b, (c*d - f*a).H]
    xy = c_[(b*d - e*a).H, (c*d - f*a).H]
    n = sqrt(wz*wz.H)
    m = sqrt(xy*xy.H)
    if n[0,0] == 0: return (Aout, Bout, Qout, Zout)
    wz = solve(n, wz)
    xy = solve(m, xy)
    wz = r_[ wz, \
            c_[-wz[:,1].H, wz[:,0].H]]
    xy = r_[ xy, \
         c_[-xy[:,1].H, xy[:,0].H]]
    Aout[ix:ix+2, :] = xy * Aout[ix:ix+2, :]
    Bout[ix:ix+2, :] = xy * Bout[ix:ix+2, :]
    Aout[:, ix:ix+2] = Aout[:, ix:ix+2] * wz
    Bout[:, ix:ix+2] = Bout[:, ix:ix+2] * wz
    Zout[:, ix:ix+2] = Zout[:, ix:ix+2] * wz
    Qout[ix:ix+2, :] = xy * Qout[ix:ix+2, :]
    return (Aout, Bout, Qout, Zout)
#
def qzdiv(stake, A2, B2, Q, Z):
    Aout = A2.copy(); Bout = B2.copy(); Qout = Q.copy(); Zout = Z.copy()
    n, jnk = A2.shape
    # remember diag returns 1d
    root = mat(abs(c_[diag(A2)[:,newaxis], diag(B2)[:,newaxis]]))
    root[:,1] /= where(root[:,0]<1e-13, -root[:,1], root[:,0])
    for i in range(1,n+1)[::-1]:        # always first i rows, decreasing
        m = None
        for j in range(1,i+1)[::-1]:    # search backwards in the first i rows
            #print root.shape
            #print n, i, j
            #print 'i,j in qzdiv', i,j
            if root[j-1,1] > stake or root[j-1,1] < -0.1:
                m = j                   # get last relevant row
                break
        if m == None: return (Aout, Bout, Qout, Zout)
        for k in range(m,i):            # from relev. row to end of first part
            (Aout, Bout, Qout, Zout) = qzswitch(k, Aout, Bout, Qout, Zout)
            root[k-1:k+1, 1] = root[k-1:k+1, 1][::-1]
    return (Aout, Bout, Qout, Zout)


def solab(A,B,n_s):
    """Port of Paul Klein's linear solver 
    (this is based on solabalt, which uses qzdiv as opposed to the Matlab ordqz).

    Purpose: Solves for the recursive representation of the stable solution to a system of linear difference equations.
    
    Inputs: Two square matrices A and B and a natural number n_s

    A and B are the coefficient matrices of the difference equation
    A*x(t+1) = B*x(t)
    where x(t) is arranged so that the state variables come first, and n_s is the number of state variables.

    Outputs: the decision rule f and the law of motion p. If we write
    x(t) = [k(t);u(t)] where k(t) contains precisely the state variables, then
    u(t)   = f*k(t) and
    k(t+1) = p*k(t)
    
    Calls: qz, qzdiv"""

    [s,t,q,z] = qz(A,B,output = 'complex')
    [s, t, q, z] = qzdiv(1., s,t,q,z)


    z21 = z[n_s:,:n_s]
    z11 = z[:n_s,:n_s]

    if np.linalg.matrix_rank(z11) < n_s:
        raise Exception("Invertibility condition violated")




    z11i = solve(z11,np.eye(n_s))
    s11 = s[:n_s,:n_s]
    t11 = t[:n_s,:n_s]

    if np.abs(t[n_s-1,n_s-1])>np.abs(s[n_s-1,n_s-1]) or np.abs(t[n_s,n_s])<np.abs(s[n_s,n_s]):
       raise Exception("Wrong number of stable eigenvalues")


    dyn = solve(s11,t11)

    f = np.dot(z21,z11i).real
    p = np.dot(z11,np.dot(dyn,z11i)).real

    return f,p


def Perturb(f1, g1):
    """Computes a Taylor approximation of decision rules, given the supplied derivatives.
    The original system is assumed to be in the the form:
    .. math::
        E_t f(s_t,x_t,s_{t+1},x_{t+1})
        s_t = g(s_{t-1},x_{t-1}, \\lambda \\epsilon_t)
    where :math:`\\lambda` is a scalar scaling down the risk.

    The solution is a triplet of matrices :math:`A, B, C` such that:
    .. math::
        s_t = A  s_{t-1} + B \\epsilon_t
        x_t = C  s_t
    The user supplies, a list of derivatives of f and g.
    :param f1: derivatives of f 
    :param g1: derivatives of g 

    """

    import numpy as np
    from numpy.linalg import solve


    n_x = f1.shape[0]           # number of controls
    n_s = f1.shape[1]/2 - n_x   # number of states
    n_e = g1.shape[1] - n_x - n_s
    n_v = n_s + n_x

    f_s = f1[:,:n_s]
    f_x = f1[:,n_s:n_s+n_x]
    f_snext = f1[:,n_v:n_v+n_s]
    f_xnext = f1[:,n_v+n_s:]

    g_s = g1[:,:n_s]
    g_x = g1[:,n_s:n_s+n_x]
    g_e = g1[:,n_v:]

    A = np.row_stack([
        np.column_stack( [ np.eye(n_s), np.zeros((n_s,n_x)) ] ),
        np.column_stack( [ -f_snext    , -f_xnext             ] )
    ])
    B = np.row_stack([
        np.column_stack( [ g_s, g_x ] ),
        np.column_stack( [ f_s, f_x ] )
    ])

    C,A = solab(A,B,n_s)
    B = g_e

    return A, B, C

## Tests ##

# ns = 2
# nx = 2
# ne = 1

# nv = ns + nx

# s = np.random.rand(ns,1)
# x = np.random.rand(nx,1)



# def f(s,x,sn,xn):
#     return np.vstack((s[0]*s[1], xn[0]*x[1]*sn[1] ))


# def g(s,x,e):
#     return np.vstack((s[0]*s[1], x[0]*x[1] * e[0] ))

# f1,g1 = Linearize(f,g,s,x,ne)


# def Test_derivs(s,x,sn,xn):

#     e = np.zeros((ne,1))

#     f1 = np.zeros((nx,2*nv))
#     f1[0,:] = np.hstack([s[1], s[0], np.zeros(nx+ns+nx)])
#     f1[1,:] = np.hstack([np.zeros(ns), 0, xn[0]*sn[1], 0, xn[0]*x[1], x[1]*sn[1], 0 ])

#     g1 = np.zeros((ns,nv+ne))
#     g1[0,:] = np.hstack([s[1], s[0], np.zeros(nx+ne)])
#     g1[1,:] = np.hstack([np.zeros(ns), x[1]*e[0],x[0]*e[0],x[0]*x[1]])
    
#     return f1, g1


# f12,g12 = Test_derivs(s,x,s,x)
# assert(np.allclose(f1.flatten(),f12.flatten()))
# assert(np.allclose(g1.flatten(),g12.flatten()))

