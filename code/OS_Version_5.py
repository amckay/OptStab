print 'Version 5 -- more volatile business cycle -- raise std of tfp and mon pol by 25%'
execfile('OS_Version_1.py')

shock_covar[0:2,0:2] *= 1.25**2
