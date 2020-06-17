import numpy as np
import os


OutDir = os.path.join(os.path.pardir,'data','results')


def ShowResults(version):

    with np.load(os.path.join(OutDir,'SearchResults_v'+str(version)+'.npz')) as X:
        OptCyc = X['OptCyc']
        OptSS = X['OptSS']




    print 'with cycles = {0:.{1}f}'.format(OptCyc[0],3)
    print 'without cycles = {0:.{1}f}'.format(OptSS[0],3)
    print ''

print('')
print('-----optimal b with zeta x 2 -----')
ShowResults(version = 4)

print('-----optimal b with wage response to benefits -----')
ShowResults(version = 2)
