from OptStab import *


w = M.wbar
Y,A,S,J,MM,h,q,Vn,u = M.SteadyState()[[iY,iA,iS,iJ,iM,ih,iq,iVn,iu]]

def g(Vn,b,Mode):
    q = (MM*Vn)**(1.0/kappa)
    Vn2 = (-log(b) -   h**(1+gamma)/(1+gamma) + psy )/(1-beta*(1-upsilon)*(1-kappa/(1+kappa)*q*MM) )
    u = upsilon*(1 - q * MM) / (q*MM + upsilon*(1 - q * MM) )
    if Mode.lower() == 'residuals':
        return np.array((Vn2 - Vn))
    elif Mode.lower() == 'values':
        return np.array((Vn, q, u))
    else:
        raise Exception('Unrecognized value for Mode')

b0 = M.b
sol = optimize.root(lambda x: g(x,b0,'residuals'),Vncalib,jac=False)
x0 = g(sol.x,b,Mode = 'values')

eps = 0.002
b1 = b + eps
sol = optimize.root(lambda x: g(x,b1,'residuals'),Vncalib,jac=False)
x1 = g(sol.x,b,Mode = 'values')

print("Micro-Elas = {0}".format((x1[2]-x0[2])/eps*b0/x0[2])    )

print("Elas of x wrt b")
print( (x1-x0)/eps*b0/x0 )

print("Deriv of x wrt b")
print( (x1-x0)/eps )


# theoretical
VnT = -log(b) / (1-beta*(1-upsilon)*(1-kappa/(1+kappa) *q*MM) )
qT = (q*MM*Vn)**(1./(kappa+1))


dudq = -MM*u**2/upsilon/(1-q*MM)**2
dVndq = log(b)*beta*(1-upsilon)*kappa/(1+kappa)*MM / (1-beta*(1-upsilon)*(1-kappa/(1+kappa)*q*MM))**2
dqdb = -1./kappa * q**(1-kappa) *MM / (1-beta*(1-upsilon)*(1-kappa/(1+kappa)*q*MM)) / b / (1-q**(1-kappa)*dVndq/kappa*MM)
#dqdb doesn't match exactly
MicroElasticity_Implied_a = b/ubar * dudq * dqdb
MicroElasticity_Implied_b = b/ubar * dudq * (x1[1]-x0[1])/eps


dVndb_a = -1.0/(1-beta*(1-upsilon)*(1-kappa/(1+kappa)*q*MM)) / b + dVndq *dqdb
dVndb_b = -1.0/(1-beta*(1-upsilon)*(1-kappa/(1+kappa)*q*MM)) / b + dVndq * (x1[1]-x0[1])/eps
