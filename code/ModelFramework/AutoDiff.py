import numpy as np
from ad.admath import exp as adexp
from ad.admath import log as adlog
from ad import adnumber, jacobian



def exp(x):
    """wrapper for ad version of exp"""
    return np.asarray(adexp(x))

def log(x):
    """wrapper for ad version of log"""
    return np.asarray(adlog(x))




def DoJac(f,X):
	"""Takes a function and a numpy array.
	Returns (J,Y) where Y = f(X) and J is the Jacobian."""

	Xad = adnumber(X)

	Yad = f(Xad)

	J = np.array(jacobian(Yad.flatten(),Xad.flatten()))

	return J, Yad.astype(float)


