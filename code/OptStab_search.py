# Command line options
# GridSearch    if supplied, do search on a grid
# NumProc N    Use N processors (default serial processing)
# CyclicalPolicy {0,1,2}:
#   0 -- (default)  bCyc = 0.0
#   2 -- bCyc = 2.1122  (From data)
#   99 -- bCyc optimized

from OptStab import *
from ModelFramework import EDS
import os
from scipy.optimize import brentq
import sys
import copy

OutDir = os.path.join(os.path.pardir,'data','results')

bCycFromData = 1.8757


WelfareFindices = iVRange-M.nX

# parse command line arguments
GridSearch = CLArgs.GridSearch
CyclicalPolicy = CLArgs.CyclicalPolicy
NumProc = CLArgs.NumProc
assert  NumProc < 29, "invalid number of processors"

print 'GridSearch = {0}'.format(GridSearch)
print 'CyclicalPolicy = {0}'.format(CyclicalPolicy)
print 'Number of processors: {0}'.format(NumProc)


def MWMSolve(M,WM,DoNonLinear = True, FindGrid = True, verbose = False):
    # initialize with pertrubation solution
    # we overload exp and log using functions that can handle ad
    # and then we switch back to the numpy versions, which are faster
    M.exp = adexp
    M.log = adlog
    PertSol = M.Perturb(SetF = True)
    M.exp = exp
    M.log = log

    # Solve the model by constructing an EDS grid
    try:
        if DoNonLinear:
            if FindGrid:
                SS = {'state0': M.SteadyState()[:M.nX][np.newaxis].T,'nrep':10}
                grid, Xsim = M.SolveEDS(T/10,nGrid,damp = [0.995,0.995],verbose = verbose, PlotGrids = False, SimStart = SS)
            else:
                M.Solve(Sgrid,damp = [0.995,0.995],verbose = verbose)

        else:
            Xsim = M.Simulate(T)
            THIN = max(0.05, 2 * float(nGrid)/float(T))
            EDS_points_tol = 2
            grid = ModelFramework.EDS.GetGrid(Xsim[:M.nX],THIN,nGrid,EDS_points_tol)[0]


        #create a copy of the coefficient matrix
        WM.F.coef = np.array(M.F.coef)
        WM.F.fitted = True

        #add the extra state
        newstates = np.setdiff1d(WelfareStates,DynamicsStates)
        if len(newstates) > 0:
            WM.F.AddStates(len(DynamicsStates),np.searchsorted(DynamicsStates,newstates) + range(len(newstates)))

            # redo the grid to include the new states
            if FindGrid:
                THIN = max(0.1, 2 * float(nGrid)/float(T))
                EDS_points_tol = 2
                grid_a = EDS.GetGrid(Xsim[WM.interpstates],THIN,nGrid,EDS_points_tol)[0].T
                noninterpstates = np.setdiff1d(range(WM.nX),WM.interpstates)
                grid = np.zeros((WM.nX,grid_a.shape[1]))
                grid[WM.interpstates] = grid_a
                grid[noninterpstates] = Xsim[noninterpstates,:grid.shape[1]]

                # add some Sobol points
                GridWidth = 1.5 * ( grid.max(axis = 1) - grid.min(axis = 1) )
                Sob = ModelFramework.Drivers.SobolGrid(grid.mean(axis = 1),GridWidth, nSobol = 20)
                grid = np.hstack((grid,Sob))

        if not FindGrid:
            grid = Sgrid


        #solve for welfare
        WM.SolvePolicyIt(grid,WelfareFindices)

        if FindGrid:
            return grid


    except RuntimeError:
        if DoNonLinear:
            warnings.warn('Non-linear solution failed for bSS = {0}, bCyc = {1}, tau = {2}.  Trying linear solution.'.format(M.bSS,M.bCyc,M.tau))
            return MWMSolve(M,WM,DoNonLinear = False, FindGrid = FindGrid)# try again without non-linear solution
        else:  # we have already tried DoNoLinear = False so now we just fail.
            raise RuntimeError('Could not solve the model for bSS = {0}, bCyc = {1}, tau = {2}.'.format(M.bSS,M.bCyc,M.tau))




def testfcn(b,tau,M,WM,WelfareFindices,DoCycles=True,DoMoments=False,DoNonLinear = True,bCyc = 0.0):


    M.Set_b_tau(b,tau,bCyc=bCyc)
    WM.Set_b_tau(b,tau,bCyc=bCyc)


    if DoMoments or DoCycles:

        try:
            MWMSolve(M,WM,DoNonLinear=DoNonLinear, FindGrid = False)
        except RuntimeError:
            if not DoMoments:
                return -100000.0
            else:
                raise RuntimeError('Trouble solving the model')

        Welfare_Cyc = WM.F(WM.SteadyState()[WM.interpstates][np.newaxis].T)[WelfareFindices].flatten()
        Mom = M.GetMoments()

    if DoMoments or not DoCycles:
        Welfare_NoCyc = WM.SteadyState()[WM.nX:WM.nXY][WelfareFindices].flatten()



    if DoMoments:
        return (Welfare_Cyc, Welfare_NoCyc, Mom['mn'], np.sqrt(np.diag(Mom['cov'])), M.Observation(M.SteadyState()))
    elif DoCycles:
        return Welfare_Cyc[3]
    else:
        return Welfare_NoCyc[3]





# -------------------------
#      Script to solve
# -------------------------

if __name__ == "__main__":


    # Create the welfare model object
    WF = PI.Poly(poly_ord)
    interpstates = WelfareStates  # use all of the states
    WM = OptStabDynamics(VarNums,shock_covar,WF,IncProc,b,tau,nquad,interpstates)


    #set up the grids for Sobol
    fname = os.path.join(OutDir,'GridBounds_v'+str(VERSION)+'.npz')
    with np.load(fname) as X:
        # gmin = X['gmin']
        # gmax = X['gmax']
        savedgrid = X['savedgrid']

    Sgrid_a = EDS.GetGrid(savedgrid[WM.interpstates],1.0,nGrid,2,ElimLowProbPoints = False)[0].T
    Sgrid = np.zeros((WM.nX,Sgrid_a.shape[1]))
    Sgrid[WM.interpstates] = Sgrid_a
    noninterpstates = np.setdiff1d(range(WM.nX),WM.interpstates)
    Sgrid[noninterpstates] = savedgrid[noninterpstates,:Sgrid.shape[1]]

    # add some Sobol points
    GridWidth = Sgrid.max(axis = 1) - Sgrid.min(axis = 1)
    Sob = ModelFramework.Drivers.SobolGrid(Sgrid.mean(axis = 1),GridWidth, nSobol = 20)
    Sgrid = np.hstack((Sgrid,Sob))


    if GridSearch:
        print 'Solving on grid'

        b_grid, tau_grid = np.meshgrid(np.arange(bmin,bmax,bstep), np.arange(taumin,taumax,taustep))
        btaulist = zip(b_grid.flatten(),tau_grid.flatten())

        def Wrap4GridSearch(btau):
            try:
                return testfcn(btau[0],btau[1],M,WM,WelfareFindices,
                    DoCycles = True, DoMoments = True, DoNonLinear = True)
            except (RuntimeError, ValueError):
                return None

        #  get the values on the mesh grid
        if NumProc > 0:
            import multiprocessing
            p = multiprocessing.Pool(NumProc)
            res = p.map(Wrap4GridSearch, btaulist)
        else:
            res = map(Wrap4GridSearch, btaulist)


        # unpack the values
        Results = [[], [], [], [], []]
        for r in res:
            if r is not None:
                for i in range(len(Results)): # unpack the results into the list of lists
                    Results[i].append(r[i])
            else:
                for i in range(len(Results)): # unpack the results into the list of lists
                    tmp = Results[i][-1]
                    if type(tmp) is np.ndarray:
                        Results[i].append(np.NaN*np.ones(tmp.shape))
                    else:
                        Results[i].append(np.NaN)


        Signature = 'WelfareGrid_v' + str(VERSION)
        if CyclicalPolicy > 0:
            Signature += "_CyclicalPolicy_" +str(CyclicalPolicy)

        np.savez(os.path.join(OutDir,Signature + '.npz'),
            WelfareCyc=np.array(Results[0]),
            WelfareNoCyc =np.array(Results[1]),
            Means = np.array(Results[2]),
            StDevs = np.array(Results[3]),
            SteadyState = np.array(Results[4]),
            b_grid = b_grid, tau_grid = tau_grid)


        print 'Done solving on grid.'

    else:
        print 'Optimizing'

        from scipy import optimize
        res1 = optimize.minimize(lambda x: -testfcn(x[0],x[1],M,WM, WelfareFindices,
                DoCycles = False, DoMoments = False, DoNonLinear = True), [0.75,0.35],method = 'Nelder-Mead' )

        print res1



        # now with cycles
        start_time = time.time()
        if VERSION == 2:
            x0 = (0.7084008288326733, 0.28529334887090951)
        else:
            x0 = res1['x']

        if CyclicalPolicy == 0:
            res2 = optimize.minimize(lambda x: -testfcn(x[0],x[1],M,WM,WelfareFindices,
                DoCycles = True, DoMoments = False, DoNonLinear = True, bCyc = 0.0), x0,method = 'Nelder-Mead', options = {'ftol': 1e-6} )
        elif CyclicalPolicy == 2:
            res2 = optimize.minimize(lambda x: -testfcn(x[0],x[1],M,WM,WelfareFindices,
                DoCycles = True, DoMoments = False, DoNonLinear = True, bCyc = bCycFromData), x0,method = 'Nelder-Mead', options = {'ftol': 1e-6} )
        elif CyclicalPolicy == 99:
            x0 = np.hstack((x0,35.0))
            res2 = optimize.minimize(lambda x: -testfcn(x[0],x[1],M,WM,WelfareFindices,
                DoCycles = True, DoMoments = False, DoNonLinear = True, bCyc = x[2]), x0,method = 'Nelder-Mead', options = {'ftol': 1e-6} )
        else:
            raise RuntimeError("CyclicalPolicy option not recognized")


        print res2

        Signature = 'SearchResults_v' + str(VERSION)
        if CyclicalPolicy > 0:
            Signature += "_CyclicalPolicy_" +str(CyclicalPolicy)

        np.savez(os.path.join(OutDir,Signature + '.npz'),OptSS = res1['x'], OptCyc = res2['x'])

        print 'Elapsed time = {}'.format(time.time()-start_time)
