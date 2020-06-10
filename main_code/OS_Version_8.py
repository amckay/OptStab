print 'Version 8 -- Flexible Prices'
execfile('OS_Version_1.py')

# this is from version 1 with the MP shocks removed
shock_covar  = np.diag( np.diag(shock_covar)[[0,2]] )




WelfareStates = range(4)  # don't include ulalg among states (too little variation under flexible prices)
