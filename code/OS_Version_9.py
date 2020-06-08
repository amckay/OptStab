print 'Version 9 -- Bigger Demand Shocks'
execfile('OS_Version_1.py')
# to get these: edit OptStab_calibrate.py to change the TechnologyVarShare to 0.25
shock_covar  = np.diag((0.04097536,    0.01453296 , 2.14196937))/1000.0
