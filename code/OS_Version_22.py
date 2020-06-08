print 'Version 22 -- Set monetary shock to 50% st dev'
execfile('OS_Version_1.py')

shock_covar[1,1] *= 0.5**2
