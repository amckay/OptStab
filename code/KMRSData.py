# Data from "Gross Worker Flows over the Business Cycle," by
# Per Krusell, Toshihiko Mukoyama, Richard Rogerson, and Aysegul Sahin.
# American Economic Review, 107 (11): 3447-3476, November 2017.

import numpy as np


# transition matrix between E, U, N  from Table 2
# (i,j) element is from i to j
EUN_Trans = np.array([[0.972, 0.014, 0.014], [0.228, 0.637, 0.135], [0.022, 0.021, 0.957]])
std_fEU = 0.089  # standard deviation of EU flow from Table 3


# get rid of flows to from N
EU_Trans = EUN_Trans[:2,:2].T
EU_Trans = EU_Trans / EU_Trans.sum(axis=0)  # condition on not going to N

# Convert to quarterly rates
EU_Trans = np.linalg.matrix_power(EU_Trans,3)

# Go back to (i,j) element is from i to j
EU_Trans = EU_Trans.T


# Alternatively
A = np.linalg.matrix_power(EUN_Trans.T,3).T
A/A[:,:2].sum(axis = 1)
# and then the job separation rate  is 3.8%
