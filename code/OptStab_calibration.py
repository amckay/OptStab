
#################
# calibration
#################
with open('../data/calibration/Moments.txt','r') as infile:
    TargetMoments = json.load(infile)

delta = 1.0/200.0
theta = 1.0/3.5
omegapi = 1.6578089766
omegau = 0.13321923864

ucalib = ubar = TargetMoments['meanu']

from KMRSData import EU_Trans, std_fEU
qM = EU_Trans[1,0]
upsilon  = ubar * qM /(1-qM - ubar * (1-qM))
fEU = upsilon*(1-qM)

tau = 0.151
b = 0.81
gamma = 2.
gamma_0 = 1
mu = 1.2
WageElasticity = 0.0
GCElas = 1.0



#GonY  = chi/(1+chi)
GonY = TargetMoments['meanGY']
chi = GonY/(1-GonY)
SSGRatio = chi
# xbar normalized to zero
# IP_u2x = np.array((0.0, 6.95, 0.0))  # see IncomeProcessMcKay2015.py
# IPP = {'sigma': (0.0143,    0.1041,0.1041), 'p': (0.,0.0526,0.0526),
#        'mu_0': (0.,.3071 ,  -0.2508), 'mu_x': (0.,1.,1.), 'neps': (5,7,7)}

if VERSION == 2:
    IncProcVersion = 2
    WageElasticity = 0.0005  # test with  python CheckElasticityofuTob.py ver 2

    # coefficeints of u_{v2}(b) - u_{v1}(b)
    #uadjustcoef = a - b, where a comes from python CheckElasticityofuTob.py ver 2 and b from python CheckElasticityofuTob.py ver 1
    uadjustcoef = np.array([  0.08560233, -0.06356474,  0.05640003])-np.array([ 0.09339892, -0.11133585,  0.08997821])

else:
    IncProcVersion = 1
    uadjustcoef = np.zeros(3)


IP_u2x = np.array((0.0, 16.73, 0.0))
IPP = {'sigma': (0.0403,0.0966,0.0966), 'p': (0.,0.00727,0.00727),
       'mu_0': (0.,0.266,-0.184), 'mu_x': (0.,1.,1.), 'neps': (5,7,7)}

IncProc = IncomeProcess(IPP,IP_u2x,ucalib,IncProcVersion,uadjustcoef)
IP_Moments = IncProc.get_Moments(ucalib,ucalib,tau,b)
beta  = 1/1.03**0.25/(((1-fEU) + fEU/b) * float(IP_Moments[2]))


Phi = 1.79
zeta_1 = 1.0
zeta_2 = (1.+gamma) / Phi
psi_1 = 0.03088
psi_2 = 1.0


rhoA = 0.9
rhoMP = 0.9
rhoG = 0.9


MicroElasticity = 0.5  # see Landais et al. NBER working paper page 44


if not VERSION == 15:
    print 'Initializing OptStab with VERSION = {0}'.format(VERSION)
    execfile('OS_Version_' + str(VERSION) + '.py')




def CalibKappa(X,Mode = 'Residuals'):
    # calibration of kappa--see notebook for 2/17/16, 2/15/18
    kappa, Y, h = X # unpack
    A = S = 1.0  # steady state values


    Vn = -log(b) / (1-beta*(1-upsilon)*(1-kappa/(1+kappa) *qM) )
    q = (qM*Vn)**(1./(kappa+1))
    M = qM/q

    dudq = -M*ubar**2/upsilon/(1-qM)**2
    dVndq = log(b)*beta*(1-upsilon)*kappa/(1+kappa)*M / (1-beta*(1-upsilon)*(1-kappa/(1+kappa)*qM))**2
    dqdb = -1./kappa * q**(1-kappa) *M / (1-beta*(1-upsilon)*(1-kappa/(1+kappa)*qM)) / b / (1-q**(1-kappa)*M*dVndq/kappa)
    MicroElasticity_Implied = b/ubar * dudq * dqdb

    w = A/mu - upsilon*psi_1*M**psi_2/h
    H = 1-ubar-(1-upsilon)*(1-ubar)
    J = psi_1 * M**psi_2 * H
    h2 = ((1-tau) * w * S  / A  * Y/ (Y-J))**(1.0/(1+gamma))
    Y2 = A * h2 * (1-ubar)

    if Mode.lower() == 'residuals':
        return np.array((MicroElasticity_Implied - MicroElasticity, Y2 - Y, h2 - h))
    elif Mode.lower() == 'values':
        return kappa,  w, Y, J, h, M, Vn
    else:
        raise Exception('Mode not recognized.')



sol = optimize.root(CalibKappa,np.array((6.0, 1.0, 1.0)),jac=False)
kappa, wcalib, Ycalib, Jcalib, hbar, Mcalib, Vncalib = CalibKappa(sol.x,Mode = 'values')
bcalib = b

if VERSION == 15:
    print 'Initializing OptStab with VERSION = {0}'.format(VERSION)
    execfile('OS_Version_' + str(VERSION) + '.py')


psy = hbar**(1+gamma)/(1+gamma) # same assumption as Baily



#recruiting cost as a share of quarterly pay.  Target = 0.03
print 'recruiting cost as a share of quarterly pay = {0}'.format((psi_1*Mcalib**psi_2)/(wcalib*hbar))
print 'calibrated value for 1/kappa = {0}'.format(1.0/kappa)
