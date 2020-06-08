print 'Version 27 -- Set rhoG to 0.8'
execfile('OS_Version_1.py')


shock_covar[2,2] *= (1-0.8**2)/(1-rhoG**2)
rhoG = 0.8
