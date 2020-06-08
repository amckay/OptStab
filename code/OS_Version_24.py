print 'Version 24 -- Set G shock to 50% st dev'
execfile('OS_Version_1.py')

shock_covar[2,2] *= 0.50**2
