import numpy as np
from ModelFramework.CompEcon import qnwnorm
from ModelFramework.Support import NaNInfError
import ad
from numpy import exp, log



class IncomeProcess:

    def __init__(self,IPP,IP_u2x,ucalib,IncProcVersion,uadjustcoef = None):
        #check income process parameters
        assert(len(IPP['sigma']) == len(IPP['neps']))
        assert(len(IPP['p']) == len(IPP['neps']))
        assert(len(IPP['mu_0']) == len(IPP['neps']))
        assert(len(IPP['mu_x']) == len(IPP['neps']))


        sigma = np.array(IPP['sigma'])[np.newaxis].transpose()

        p = np.array(IPP['p'])
        p[0] = 1.0 - p[1:].sum()
        p = p[np.newaxis].transpose()

        # mu_0 = np.array((0.,0.3550,-0.2989))[np.newaxis].transpose()
        # mu_x = np.array((0.,1.,1.))[np.newaxis].transpose()

        mu_0 = np.array(IPP['mu_0'])[np.newaxis].transpose()
        mu_x = np.array(IPP['mu_x'])[np.newaxis].transpose()



        # QUADRATURE OVER INCOME SHOCKS
        neps = IPP['neps']
        cumsumneps = np.cumsum(np.hstack((0,neps)))
        epsgrid = np.array([])
        epsprob = np.array([])
        for i in range(len(neps)):
            eg,ep = qnwnorm(neps[i],0.,sigma[i]**2)
            epsgrid = np.hstack((epsgrid,eg))
            epsprob = np.hstack((epsprob,p[i]*ep))

        epsgrid = epsgrid[np.newaxis].transpose()
        epsprob = epsprob[np.newaxis].transpose()


        self.neps = neps
        self.p = p
        self.mu_0 = mu_0
        self.mu_x = mu_x
        self.epsgrid = epsgrid
        self.epsprob = epsprob
        self.cumsumneps = cumsumneps
        self.sigma = sigma

        self.IP_u2x = IP_u2x
        self.ucalib = ucalib
        self.IncProcVersion = IncProcVersion
        self.uadjustcoef = uadjustcoef
        assert IncProcVersion == 1 or IncProcVersion == 2

    def __ShockGrid__(self, x ):
        """Returns the quadrature grid and weights for a given level of shock
        the grid is for log(epsilon)

        epsgrid, epsprob = ShockGrid( x )

        Inputs:
        x -- scalar level of risk  (can be a row vector)



        Outputs:
        epsgrid -- grid for eps realizations (dimension sum(neps) by len(x))
        epsprob -- probability weights
        """


        MU  = self.__get_MU__(x )

        epsgrid1 = MU + self.epsgrid

        return epsgrid1, self.epsprob




    def __get_MU__(self, x ):
        """Returns the three central moments of the mixture components for the distribtion of eps."""

        MU = self.mu_0 - self.mu_x * x

        # mubar = -log((self.p * exp(MU + self.sigma**2/2.0)).sum(axis=0))
        # MU = mubar + MU


        # format MU by repeating the rows the appropriate number of times
        if x.__class__ == np.ndarray :
            width = len(x)
            isad = x[0].__class__ == ad.ADV
        else:
            width = 1
            isad = x.__class__ == ad.ADV

        # if we are doing automatic differentiation we need to set the data type to object
        # if we do not do this, then assigning into MU_formatted[a:b,:] will cast the
        # ad objects to floats
        if isad:
            MU_formatted = np.zeros((np.sum(self.neps), width),dtype = object)
        else:
            MU_formatted = np.zeros((np.sum(self.neps), width))



        for i in range(len(self.neps)):
            a = self.cumsumneps[i]
            b = self.cumsumneps[i+1]
            MU_formatted[a:b,:] = MU[i,:]


        mubar = -log( (self.epsprob * exp(MU_formatted + self.epsgrid)).sum(axis = 0) )
        MU_formatted = MU_formatted + mubar

        if any( np.isnan(MU_formatted.reshape(-1)) | np.isinf(MU_formatted.reshape(-1))):
            raise NaNInfError("Bad values for MU_formatted.")

        return MU_formatted


    def get_Moments(self, u, ubar, tau, b = None):
        """Returns the moments of interest of the epsilon distribution.

        In:
        u   (float or row vector)  unemployment rate
        tau  float

        Out: M1, M2, M3  (each of same shape as x)
        M1 is E [ epsilon**(1-tau)]
        M2 is E [ log epsilon ]
        M3 is E [ epsilon**(tau -1)]
        M4 is E [ epsilon**(1-tau) * log epsilon]
        """

        # Two options here: 1) x= u - ucalib  2) x = (1-u)/(1-ubar)
        if self.IncProcVersion == 1:
            x = np.polyval(self.IP_u2x, u - self.ucalib)
        else:
            assert b is not None
            x = np.polyval(self.IP_u2x, u - self.ucalib - np.polyval(self.uadjustcoef,b))

        g, p = self.__ShockGrid__( x )

        M1 = np.dot(p.flatten(), exp( (1.0-tau) * g ))
        M2 = np.dot(p.flatten(), g )
        M3 = np.dot(p.flatten(), exp( (tau-1.0) * g ))
        M4 = np.dot(p.flatten(), exp( (1.0-tau) * g ) * g)

        return M1, M2, M3, M4
