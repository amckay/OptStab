print 'This is a version with 5% of GDP more G in steady state and acyclical G.'
execfile('OS_Version_1.py')


print('GonY = {0}'.format(GonY))
print('chi = {0}'.format(chi))
target_GY_ratio =  GonY + 0.05
SSGRatio = target_GY_ratio/(1-target_GY_ratio)
print('SSGRatio = {0}'.format(SSGRatio))
