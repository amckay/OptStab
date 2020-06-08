print 'Version 25 -- Set rhoA to 0.8'
execfile('OS_Version_1.py')


shock_covar[0,0] *= (1-0.8**2)/(1-rhoA**2)
rhoA = 0.8
