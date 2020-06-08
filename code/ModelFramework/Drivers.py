from __future__ import print_function
import numpy as np
import EDS
import sobol
from CompEcon import qnwnorm
import DetSimul as DS
from math import ceil
import sys
import warnings

class Model(object):
    """A base class for macroeconomic models.

    We use some conventions:
        - X[:nX] = State variables
        - X[nX:(nX+nY)] = Forward-looking (approximated) variables
        - X[(nX+nY):] = Static variables (solved exactly from X and Y)
        - Functions all accept and return numpy arrays.
        - X is column vector or array with columns referring either to grid points or points in time.



    The child class must implement four functions:
        - SteadyState -- returns the steady state for X
        - StateTrans -- Given values of X, and eps'  it returns states in X'
        - ForwardLooking -- Given X,X' it returns the values for X[nX:(nX+nY)] implied by the forward looking equations.
        - Static -- Given X[:nX+nY] it returns X[nX+nY:]

    Once the model has been specified there are several ways of solving it:
        - Perturb does a first order linear approximation and solves using Klein's method.
        - Solve does a global solution for a grid supplied by the user.
        - SolveEDS does a global solution by first creating an EDS grid.
        - DetSimul does a perfect foresight transition for a given set of exogenous shocks.
    """

    def __init__(self,VarNums,shock_covar,F,nquad,interpstates = None):

        # num state, forward-looking, static
        self.nX, self.nY, self.nZ  = VarNums
        if interpstates == None:
            self.interpstates = range(self.nX)
        else:
            self.interpstates = interpstates

        self.nXY = self.nX + self.nY
        self.nXYZ = sum(VarNums)

        # shock covariance matrix
        self.shock_covar = shock_covar

        #initial guess using steady state values
        Xstar = self.SteadyState()

        # interpolation object
        self.F = F
        self.F.FitXY(Xstar[self.interpstates][np.newaxis].T,Xstar[self.nX:self.nXY][np.newaxis].T)

        # get quadrature grid
        quad_x, quad_w = qnwnorm(nquad,0.0,shock_covar)
        if quad_x.ndim == 1:
            quad_x = quad_x[np.newaxis].T
        self.quad = { 'n':nquad,   'x':quad_x,   'w':quad_w }

    def CheckSteadyState(self,X):
        """Returns the residuals in the steady state relationships."""
        RX = X[:self.nX] - self.StateTrans(X,np.zeros(self.shock_covar.shape[0]))
        RY = X[self.nX:self.nXY] - self.ForwardLooking(X,X)
        RZ = X[self.nXY:] - self.Static(X)

        return np.vstack((RX, RY, RZ))

    def SolveSteadyState(self,X0):
        """Solves for the steady state from an initial guess X0"""

        from AutoDiff import DoJac
        from scipy import optimize


        def f(x):
            return self.CheckSteadyState(x)

        def jac(x):
            return DoJac(f,x)[0]


        sol = optimize.root(f, X0, jac=jac, method='hybr')
        return sol.x


    def TestSteadyState(self):
        """This is a consistency check on the model equations and the steady state function"""


        SSR = self.CheckSteadyState(self.SteadyState()[np.newaxis].T)
        if max(abs(SSR)) > 1e-6:
            for i in range(SSR.shape[0]):
                print('{0} --- {1}'.format(i,SSR[i]))

            raise Exception('steady state residuals are too large')


    def ForwardLookingExpectation(self,State):
        """Takes an expectation over future shocks"""

        assert(State.ndim == 2)
        assert(State.shape[0] == self.nX)

        X = np.zeros((self.nXYZ,State.shape[1]))
        X[:self.nX] = State
        X[self.nX:self.nXY] = self.F(X[self.interpstates])
        X[self.nXY:] = self.Static(X)

        Xprime = np.zeros(X.shape)

        Yimplied = 0.0

        for quad_i in range(np.prod(self.quad['n'])):
            Xprime[:self.nX] = self.StateTrans(X,self.quad['x'][quad_i])
            Xprime[self.nX:self.nXY] = self.F(Xprime[self.interpstates])
            Xprime[self.nXY:] = self.Static(Xprime)

            Yimplied += self.quad['w'][quad_i] * self.ForwardLooking(X,Xprime)

        return  Yimplied

    def Solve(self,grid,tol = 1e-5, maxit = 5000, damp = [0.96, 0.98],verbose = True,IncludeInFcheck = None):
        """Solves the model by iterating on the forward looking equations.
        In:
        grid   d x n   arrary of grid points
        tol   float  convergence criterion for sum of absolute relative changes
        maxit   int   max number of iterations
        damp   list of two floats  in [0,1]
            damp[0] is initial damping coefficient (weight on previous solution)
            damp[1] rate of decline of damping
            damping = max(  damp[0] * damp[1]**i, 0.5)

        IncludeInFcheck  indices of forward-looking variables to include in convergence check (default is all)

        Sets:
        self.F  solution

        Uses
        self.ForwardLooking  function that gives implied values for the forward looking variables on the grid.
        self.F  function with approximate solution for the forward looking variables.

        """

        # Test for multi-colinearity problem in interpolation
        # This warning indicates a possible multi-colinearity problem in the interpolation problem.
        # This will occur if the model is specified with state variables that have little or no variance.
        cn = self.F.GetConditionNumber(grid[self.interpstates])
        if cn > 1e8:
            warnings.warn('Large condition number in Drivers.Solve. Possible multi-colinearity problem in interpolation.')



        # apply default behavior for convergence check
        if IncludeInFcheck == None:
            IncludeInFcheck = range(self.nY)

        tstval = self.F(grid[self.interpstates])[IncludeInFcheck]
        for it in range(maxit):
            dampi = max(damp[0] * damp[1]**it,0.5)
            self.F.FitXY(grid[self.interpstates],self.ForwardLookingExpectation(grid),dampi)
            tstval2 = self.F(grid[self.interpstates])[IncludeInFcheck]
            tst = (np.abs(tstval2-tstval)/(np.abs(tstval)+tol )).sum()
            if tst < 0:
                print(tstval2)
                print(tstval)
                raise RuntimeError('negative test value')

            if verbose:
                #print tst
                print(tst, end='\r')
                sys.stdout.flush()
            if tst < tol:
                break

            tstval = tstval2

        if it == maxit-1:
            print('Solve may not have converged.')



    def Simulate(self,T,state0=None,SetSeed=True):
        """Simulates the solved model

        In:
        T  number of periods
        state0  initial state vector

        Out:
        Xsim   N x T  vector of simulated states

        uses:
        self.F  solution for forward looking variables
        self.shock_covar  variance of the exogenous shock
        """




        Xsim = np.zeros((self.nXYZ,T))
        if state0.__class__ == None.__class__:
            Xsim[:,0] = self.SteadyState()
        else:
            Xsim[:self.nX,[0]] = state0
            Xsim[self.nX:self.nXY,[0]] = self.F(Xsim[[self.interpstates],[0]].T)
            Xsim[self.nXY:,[0]] = self.Static(Xsim[:,[0]])


        if SetSeed:
            np.random.seed(9028342)

        if len(self.shock_covar.shape) == 0 or self.shock_covar.shape[0] == 1:
            nshock = 1
            eps = np.sqrt(self.shock_covar) * np.random.randn(nshock,T)
        else:
            nshock = self.shock_covar.shape[0]
            eps = np.dot(np.linalg.cholesky(self.shock_covar),np.random.randn(nshock,T))


        for t in range(1,T):
            Xsim[:self.nX,[t]] = self.StateTrans(Xsim[:,[t-1]],eps[:,[t]])
            Xsim[self.nX:self.nXY,[t]] = self.F(Xsim[[self.interpstates],[t]].T)
            Xsim[self.nXY:,[t]] = self.Static(Xsim[:,[t]])

        return Xsim

    def EDSSimulate(self,T,SimStart=None):
        '''Wrapper for simulate that allows for multiple repetitions'''
        if SimStart.__class__ == None.__class__:
            Xsim = self.Simulate(T)
        else:
            Xsim = np.zeros((self.nXYZ,T*SimStart['nrep']))
            for simit in range(SimStart['nrep']):
                Xsim[:,simit*T+np.arange(T)] = self.Simulate(T,state0 = SimStart['state0'],SetSeed = simit==0)

        return Xsim

    def SolveEDS(self,T,nGrid,tol = 1e-5, maxit = 5000, damp = [0.96, 0.98],
                PlotGrids = False,verbose = True,IncludeInFcheck = None,
                SimStart = None):
        """Solves the model by iterating on the forward looking equations.
        Constructs the EDS grid to do so.
        In:
        T  number of periods
        nGrid  number of grid points in the EDS grid
        state0  initial state vector

        tol   float  convergence criterion for sum of absolute relative changes
        maxit   int   max number of iterations
        damp   list of two floats  in [0,1]
            damp[0] is initial damping coefficient (weight on previous solution)
            damp[1] rate of decline of damping
            damping = max(  damp[0] * damp[1]**i, 0.5)
        PlotGrids Boolean (default False).  If true, the new and old grids are plotted each iteration.
        SimStart   if = None then the simulation starts at the steady state
                   otherwise, SimStart dictionary
                   SimStart['state0'] is a vector of state variables from which
                   to start the simulation.
                   SimStart['nrep'] is an integer number of times to repeat the simulation.


        Out: grid, Xsim
        grid is a d by nGrid matrix
        Xsim simulated data: nvar by T matrix

        Uses:
        self.F   initial guess of function with approximate solution for the forward looking variables.
        self.ForwardLooking  function that gives implied values for the forward looking variables on the grid.
        self.StateTrans
        self.Static
        self.shock_covar  variance of the exogenous shock
        """

        EDS_MAXIT = 20
        THIN = max(0.1, 2 * float(nGrid)/float(T))
        EDS_points_tol = 2

        noninterpstates = np.setdiff1d(range(self.nX),self.interpstates)


        Xsim = self.EDSSimulate(T,SimStart)

        if any(np.isnan(Xsim.ravel())):
            raise RuntimeError('SolveEDS encountered NaN in simulation.')

        # grid, eps, PCGrid = EDS.GetGrid(Xsim[:self.nX],THIN,nGrid,EDS_points_tol)
        grid_a, eps, PCGrid = EDS.GetGrid(Xsim[self.interpstates],THIN,nGrid,EDS_points_tol)
        grid_a = grid_a.T
        grid = np.zeros((self.nX,grid_a.shape[1]))
        grid[self.interpstates] = grid_a
        grid[noninterpstates] = Xsim[noninterpstates,:grid.shape[1]]
        X = grid

        # add some Sobol points
        GridWidth = 1.5 * ( X.max(axis = 1) - X.min(axis = 1) )
        Sob = SobolGrid(X.mean(axis = 1),GridWidth, nSobol = 20)
        X = np.hstack((X,Sob))

        for eds_it in range(EDS_MAXIT):
            self.Solve(X,tol, maxit , damp ,verbose,IncludeInFcheck)
            Xsim = self.EDSSimulate(T,SimStart)
            if any(np.isnan(Xsim.ravel())):
                raise RuntimeError('SolveEDS encountered NaN in simulation.')

            # grid, neweps, newPCGrid = EDS.GetGrid(Xsim[:self.nX],THIN,nGrid,EDS_points_tol)
            grid_a, neweps, newPCGrid = EDS.GetGrid(Xsim[self.interpstates],THIN,nGrid,EDS_points_tol)
            grid_a = grid_a.T
            grid = np.zeros((self.nX,grid_a.shape[1]))
            grid[self.interpstates] = grid_a
            grid[noninterpstates] = Xsim[noninterpstates,:grid.shape[1]]


            # add some Sobol points
            GridWidth = 1.5 * ( grid.max(axis = 1) - grid.min(axis = 1) )
            Sob = SobolGrid(grid.mean(axis = 1),GridWidth, nSobol = 20)
            grid = np.hstack((grid,Sob))

            #-----check convergence-----
            # first check that mean has converged
            if verbose:
                print('check means:')
                print(np.abs(np.mean(X,axis = 1) - np.mean(grid,axis=1))/(np.std(X,axis=1) + 0.05*np.abs(np.mean(grid,axis=1))))

            converged = all(np.abs(np.mean(X,axis = 1) - np.mean(grid,axis=1))/(np.std(X,axis=1) + 0.05*np.abs(np.mean(grid,axis=1))) < 0.1)


            if converged:
                if verbose:
                    print('mean converged')

                #next check that coverage has converged in principal component space
                for i in range(newPCGrid.shape[0]):
                    mindist = np.min(EDS.Dist(PCGrid,newPCGrid[i,:]))
                    if mindist > 2*np.min((neweps,eps)):
                        # print [mindist,  2*np.min((neweps,eps))]
                        converged = False
                        break



            if converged:
                if verbose:
                    print('Grid has converged.')

                return grid, Xsim
            else:
                if verbose:
                    print('updating grid')
                    print('neweps = ' + str(neweps))


                X = grid
                eps = neweps
                PCGrid = newPCGrid


        raise RuntimeError('SolveEDS did not converge.')


    def Perturb(self, SetF = True):
        """Solves for the decision rules using a first-order perturbation.

        Returns
        The solution is a triplet of matrices :math:`A, B, C` such that:
        .. math::

            X_t &= A  X_{t-1} + B \epsilon_t \\
            Y_t &= C  X_t \\

        where :math:`X_t` are state variables :math:`\epsilon_t` are shocks, and :math:`Y_t` are forward-looking and static variables (stacked).


        If SetF = True:
        The solution is used as the object policy rule.
        That is self.F is such that :math:`Y_t = F ( X_t )`.
        """

        X = self.SteadyState()

        ns = self.nX
        nx = self.nY + self.nZ

        ne = self.shock_covar.shape[0]


        def f(s,x,sn,xn):
            return np.vstack(( x[:self.nY] - self.ForwardLooking(np.vstack((s,x)),np.vstack((sn,xn))),
                    x[self.nY:] - self.Static(np.vstack((s,x))  )  ))

        def g(s,x,e):
            return self.StateTrans(np.vstack((s,x)),e)

        import Perturbation

        x = X[self.nX:][np.newaxis].T
        s = X[:self.nX][np.newaxis].T
        f1, g1 = Perturbation.Linearize(f,g,s,x,ne)

        A, B, C = Perturbation.Perturb(f1, g1)

        if SetF:
            C2 = C[:,self.interpstates]
            X0 = SobolGrid(X[self.interpstates], np.ones(len(self.interpstates)))
            Y0 = X[self.nX:self.nXY][np.newaxis].T + np.dot(C2[:self.nY], X0 - X[self.interpstates][np.newaxis].T)
            X0 = SobolGrid(X[self.interpstates], np.ones(len(self.interpstates)))
            self.F.FitXYLinear(X0,Y0)

        return A, B, C


    def Observation(self,X):
        return X


    def GetMoments(self, T = 2000):
        """Compute endogenous moments through simulation.

        Mom = GetMoments(T)

        Observation is a function that takes three numpy arrays (X,Y,Z) and returns the observed variables.
        X is of shape num_variables_X by T and similarly for Y and Z. The returned observed variables
        should be of shape num_variables_observed by T.

        T is number of periods to simulate (default = 2000).

        Mom['mn'] is an array of means
        Mom['cov'] is the covariance matrix
        Mom['autocorr'] is an array of autocorrelations

        For all aspects of Mom, the order of variables is determined by the Observation equation.
        """

        Y = self.Observation(self.Simulate(T))

        Mom = {}
        Mom['mn'] = Y.mean(axis = 1)
        Mom['cov'] = np.cov(Y)
        Mom['autocorr'] = np.diag(np.corrcoef(Y[:,:-1],Y[:,1:]),Y.shape[0])

        return Mom


    def DetSimul(self,X0,eps,maxit = 100,Dampening = 1, DampeningThresh = 0):
        return DS.DetSimul(X0,self.StateTrans,self.ForwardLooking,self.Static, eps, self.nX,self.nY,self.nZ, maxit, Dampening, DampeningThresh )




def SobolGrid(state0,GridWidth, nSobol = 100):
    """Returns a Sobol grid centered on state0 with total width in each dimension
    given by GridWidth.  Both state0 and GridWidth should be row vectors"""

    SOBOLSKIP = 2

    # ---- generate Sobol points ----
    X = sobol.i4_sobol_generate ( len(state0), nSobol, SOBOLSKIP )

    # scale accoring to the bounds in config file
    for i in range(len(state0)):
        X[i] = state0[i] + GridWidth[i] * ( X[i] - 0.5)

    return X
