print 'Version 26 -- Set rhoMP to 0.8'
execfile('OS_Version_1.py')

shock_covar[1,1] *= (1-0.8**2)/(1-rhoMP**2)
rhoMP = 0.8
