from OptStab import *
from OptStab_search import MWMSolve
from ModelFramework import EDS
from ModelFramework.sobol import i4_sobol_generate as sobol_generate
import os

OutDir = os.path.join(os.path.pardir,'data','results')



def UpdateSavedGridBounds(grid):


    if os.path.isfile(fname):
        print 'reading from ' + fname
        with np.load(fname) as X:
            # gmin = X['gmin']
            # gmax = X['gmax']
            savedgrid = X['savedgrid']

        # gmin = np.minimum(np.min(grid,1),gmin)
        # gmax = np.maximum(np.max(grid,1),gmax)
        savedgrid = np.hstack((savedgrid,grid))
    else:
        # gmin = np.min(grid,1)
        # gmax = np.max(grid,1)
        savedgrid = grid

    print 'writing to ' + fname
    np.savez(fname,
            # gmin = gmin,
            # gmax = gmax,
            savedgrid = savedgrid  )



# -------------------------
#      Script to solve
# -------------------------

if __name__ == "__main__":


    # clean up and do record keeping
    fname = os.path.join(OutDir,'GridBounds_v'+str(VERSION)+'.npz')
    if os.path.isfile(fname):
        os.remove(fname)

    # Create the welfare model object
    WF = PI.Poly(poly_ord)
    interpstates = WelfareStates  # use all of the states
    WM = OptStabDynamics(VarNums,shock_covar,WF,IncProc,b,tau,nquad,interpstates)
    WelfareFindices = np.array((iV_a,iV_b,iV_c,iV))-WM.nX

    b_grid, tau_grid = np.meshgrid(np.array((bmin,bmax)), np.array((taumin,taumax)))


    # ---- generate Sobol points ----

    nSobolPoints = 12
    b_grid, tau_grid = sobol_generate ( m = 2, n = nSobolPoints , skip = 2 )
    b_grid =  bmin + (bmax - bmin) * b_grid
    tau_grid =  taumin + (taumax - taumin) * tau_grid

    for b,tau in np.vstack((b_grid,tau_grid)).T:

        print 'finding grid for b = {0}, tau = {1}'.format(b,tau)
        M.Set_b_tau(b,tau)
        WM.Set_b_tau(b,tau)
        grid = MWMSolve(M,WM,DoNonLinear = True, FindGrid = True)


        print 'max'
        print np.max(grid,1)
        print 'min'
        print np.min(grid,1)

        UpdateSavedGridBounds(grid)
