print 'Version 23 -- Set productivity shock to 50% st dev'
execfile('OS_Version_1.py')

shock_covar[0,0] *= 0.5**2
