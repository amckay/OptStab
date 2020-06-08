"""Functions adapted from Lilia Maliar and Serguei Maliar (MM)

Adapted by Alisdair McKay
Saturday, July 26, 2014

MM license:

Lilia Maliar and Serguei Maliar agree to make the Software for generalized 
stochastic-simulation algorithm (GSSA) accompanying the article "Numerically 
Stable and Accurate Stochastic Simulation Approaches for Solving Dynamic 
Economic Models" by Kenneth L. Judd, Lilia Maliar and Serguei Maliar, 
published in Quantitative Economics (2011), 2/2, 173-210, available to you 
under the following conditions.

1. Any work that contains some results derived from the Software or a
   modified version must: 

   a. provide a prominent description of the use of the Software in the text; 

   b. cite the article "Numerically Stable and Accurate Stochastic Simulation 
   Approaches for Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia 
   Maliar and Serguei Maliar (2011), Quantitative Economics 2/2, 173-210. 

2. The Software and any modified version may be redistributed under the 
   terms of the present license agreement. Redistributions must include 
   the original Software, this license agreement and the source code of 
   the modified version with a prominent notice stating that the Software 
   was changed and the date of the change.

3. The Software comes as is and is used at your own risk. Any modification 
   of the Software is your own responsibility.  The Authors assume no 
   obligation to provide assistance with using the Software, make no 
   warranties concerning the performance of the Software and are not 
   liable to you or any entity for any damage related to the use of the 
   Software. 

4. The use of the Software is restricted to noncommercial research and 
   educational purposes. If the Software or a modified version leads to the 
   development of materials (including software, teaching materials, patents), 
   which may be used for commercial purposes, a specific agreement must be 
   negotiated between the authors and the beneficiary. 

5. The Software is protected by copyright and other applicable laws. 


"""






import numpy as np



def Monomials_1(N,vcv):
    """Monomials_1 is a routine that constructs integration nodes and weights 
    under N-dimensional monomial (non-product) integration rule with 2N nodes; 
    see the Supplement to "Numerically Stable and Accurate Stochastic Simulation 
    Approaches for Solving Dynamic Economic Models" by Kenneth L. Judd, Lilia 
    Maliar and Serguei Maliar, (2011), Quantitative Economics 2/2, 173-210 
    (henceforth, JMM, 2011).

    This version: July 14, 2011. First version: August 27, 2009.
    -------------------------------------------------------------------------

    (n_nodes,epsi_nodes,weight_nodes) = Monomials_1(N,vcv)

    Inputs:  "N" is the number of random variables; N>=1;
    "vcv" is the variance-covariance matrix; N-by-N

    Outputs: "n_nodes" is the total number of integration nodes; 2*N;
         "epsi_nodes" are the integration nodes; n_nodes-by-N;
         "weight_nodes" are the integration weights; n_nodes-by-1
    -------------------------------------------------------------------------"""

    
    n_nodes   = 2*N       #% Total number of integration nodes

    # 1. N-dimensional integration nodes for N uncorrelated random variables with 
    # zero mean and unit variance
    # ---------------------------------------------------------------------------   
    z1 = np.zeros((n_nodes,N)) # A supplementary matrix for integration nodes; n_nodes-by-N
                       
    for i in range(N):            # In each node, random variable i takes value either
                           # 1 or -1, and all other variables take value 0
        z1[2*i:2*(i+1),i] = (1,-1)  
    # For example, for N = 2, z1 = [1 0; -1 0; 0 1; 0 -1]

    # z = z1*sqrt(N);      % Integration nodes  

    # 2. N-dimensional integration nodes and weights for N correlated random 
    # variables with zero mean and variance-covaraince matrix vcv 
    # ----------------------------------------------------------------------  
    sqrt_vcv = np.linalg.cholesky(vcv).transpose()  # Cholesky decomposition of the variance-covariance matrix
                                 
    R = np.sqrt(N)*sqrt_vcv  # Variable R; see condition (B.7) in the Supplement to JMM (2011)
                                 
    epsi_nodes = np.dot(z1,R)     # Integration nodes; see condition ((B.7) in the Supplement% to JMM (2011); n_nodes-by-N

    # 3. Integration weights
    #-----------------------
    weight_nodes = np.ones((n_nodes,1))/float(n_nodes)
    # Integration weights are equal for all integration 
    # nodes; n_nodes-by-1; the weights are the same for 
    # the cases of correlated and uncorrelated random 
    # variables


    return (n_nodes,epsi_nodes,weight_nodes)
