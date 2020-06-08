print 'Version 7 -- Smaller Demand Shocks'
execfile('OS_Version_1.py')
# to get these, edit OptStab_calibrate.py to set TechnologyVarShare to 0.75
shock_covar  = np.diag((0.12292609,   0.00484432 , 2.14196937))/1000.0
