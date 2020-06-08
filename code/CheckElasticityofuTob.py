from OptStab import *

def ssu(b,tau):
    M.Set_b_tau(b,tau)
    return M.SteadyState()[iu]

eps = 0.01
u0 = ssu(bcalib,tau)
u1 = ssu(bcalib+eps,tau)
elas = bcalib/u0 *(u1-u0)/eps

print('Steady state elasticity of u to b = {0}'.format(elas))


b_grid = np.linspace(0.70,0.9,21)
u =  np.array(map(lambda i: ssu(i,tau), b_grid  ) )
pcoef = np.polyfit(b_grid, u, 2)

print('Quadratic fits u on b:')
print(pcoef)
