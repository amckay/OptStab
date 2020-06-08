""" Polynomial interpolation functions"""

import numpy as np


def Ord_Polynomial_N(z,D):
    """  Ord_Polynomial_N constructs the basis functions of
    complete ordinary polynomial of the degrees one to four.
    -------------------------------------------------------------------------
    Inputs: "z" is the data points on which the polynomial basis functions must be constructed; dimen x n_pts;
    "D" is the degree of the polynomial whose basis functions mustbe constructed; (can be 1,...,4)

    Output:  "basis_fs" is the matrix of basis functions of a complete polynomial of the given degree """


    assert(len(z.shape) == 2)


    dimen, n_pts = z.shape # Compute the number of rows, n_pts, and the
                          # number of variables (columns), dimen, in the
                          # data z on which the polynomial basis functions
                          # must be constructed


    # initialization
    numbasis = 1
    tmp_1 = 1
    tmp_2 = 0
    for i in range(1,D+1):
        tmp_1 *= (D+tmp_2)/i
        numbasis += tmp_1

    numbasis = 1
    tmp_1 = 1
    for i in range(1,D+1):
        tmp_1 = tmp_1 * (dimen+(i-1))/i
        numbasis += tmp_1

    basis_fs = np.ones((numbasis,n_pts))


    # 1. The matrix of the basis functions of the first-degree polynomial
    # (the default option)
    # -------------------------------------------------------------------
    basis_fs[1:dimen+1] = z
    # The matrix includes a column of ones

    i = dimen   # Index number of a polynomial basis function; the first
                 # basis function (equal to one) and linear basis functions
                 # are indexed from 1 to dimen, and subsequent polynomial
                 # basis functions will be indexed from dimen+1 and on

    # 2. The matrix of the basis functions of the second-degree polynomial
    # --------------------------------------------------------------------
    if D == 2:
        for j1 in range(dimen):
            for j2 in range(j1,dimen):
                i = i+1
                basis_fs[i] = z[j1]*z[j2]

    # 3. The matrix of the basis functions of the third-degree polynomial
    # -------------------------------------------------------------------
    elif D == 3:
        for j1 in range(dimen):
            for j2 in range(j1,dimen):
                i = i+1
                basis_fs[i] = z[j1]*z[j2]
                for j3 in range(j2,dimen):
                    i = i+1
                    basis_fs[i] = z[j1]*z[j2]*z[j3]

    # 4. The matrix of the basis functions of the fourth-degree polynomial
    # -------------------------------------------------------------------
    elif D == 4:

        for j1 in range(dimen):
            for j2 in range(j1,dimen):
                i = i+1;
                basis_fs[i] = z[j1]*z[j2]
                for j3 in range(j2,dimen):
                    i = i+1
                    basis_fs[i] = z[j1]*z[j2]*z[j3]
                    for j4 in range(j3,dimen):
                        i = i+1
                        basis_fs[i] = z[j1]*z[j2]*z[j3]*z[j4]




    return basis_fs

def  PolyIntBasis(X,ord):
    """Wraps basis function call for easy changes"""

    return Ord_Polynomial_N(X,ord)

class Poly:

    def __init__(self,ord,fname = None):
        """Sets up a polynomial interpolation structure.
        s = AMpolyint_step1( X, ord,s0 )
        In:
        ord       order (1, 2, 3)
        fname     optional filename with numpy array with coefficients
        """

        self.ord = ord
        if fname is not None:
            self.coef = np.load(fname)
            self.Fitted = True
        else:
            self.Fitted = False


    def save(self,fname):
        """Save coefficients to a file"""
        np.save(fname,self.coef)

    def FitXY(self,X,Y,damping = 0.0):
        """
        Finds coefficients to fit y
        In: X    d x n  array of points on which we will fit the polynomial (ie the grid)
            Y is a matrix of functions evaluated on the grid [columns are separate functions]
            """

        if len(X.shape) == 1:
            X = X[np.newaxis]

        if len(Y.shape) == 1:
            Y = Y[np.newaxis]

        newcoef = np.linalg.lstsq(PolyIntBasis(X,self.ord).transpose(),Y.transpose(),rcond=None)[0].transpose()

        if damping > 0.0:
            self.coef = damping * self.coef + (1-damping) * newcoef
        else:
            self.coef = newcoef
        self.Fitted = True

    def FitXYLinear(self,X,Y):
        """
        Finds coefficients to fit Y setting non-linear coefficients to zero.
        In: X    d x n  array of points on which we will fit the polynomial (ie the grid)
            Y is a matrix of functions evaluated on the grid [columns are separate functions]
            """
        CoefLinear = np.linalg.lstsq(PolyIntBasis(X,1).transpose(),Y.transpose(),rcond=None)[0].transpose()
        self.coef = np.hstack((CoefLinear,np.zeros((len(Y),len(PolyIntBasis(X,self.ord))-CoefLinear.shape[1]))))
        self.Fitted = True


    def GetConditionNumber(self,X):
        """
        Finds coefficients to fit y
        In: X    d x n  array of points on which we will fit the polynomial (ie the grid)
        Out: the condition number of the basis matrix
            """
        B = PolyIntBasis(X,self.ord).transpose()
        return np.linalg.cond(B)





    def __call__( self, X ):
        """y = s(X)
        Interpolates functions at X
        In:  X   d x n  array of points to interpolate at
        Out: Y   nfunc x n array of interpolated values"""


        if len(X.shape) == 1:
            X = X[np.newaxis]

        assert(self.Fitted)

        Y = np.dot(self.coef,PolyIntBasis(X,self.ord))


        return Y



    def AddStates(self, numexistingstates, newstateindices):
        '''Reformats coef matrix to allow for new states in the state vector.
        The coefficients on the new states are set to zero.

        numexistingstates = number of states in the original state vector
        newstateindices  = list of indices of the new states (as indices of the new state vector)

        The new state vector has length n = numexistingstates + len(newstateindices)

        newstateindices must be between 0 and n without duplicates.'''


        n = numexistingstates + len(newstateindices)

        X = np.ones((n,1))
        X[newstateindices] = 0

        Z = PolyIntBasis(X,self.ord).flatten().astype(bool)
        newcoef = np.zeros((self.coef.shape[0],len(Z)))
        newcoef[:,Z] = self.coef

        self.coef = newcoef
