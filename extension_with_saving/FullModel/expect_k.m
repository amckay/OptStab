% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% Compute j-th central moment of the cross-sectional distribution of k,
% Inputs:
%   pvec:  vector of probabilities
%          with 2 columns: first column vector of k-nodes, second column vector of probabilit
%   j:     which_ moment
function EV = expect_k(pvec,j,k)
  global Params;
  if(nargin==1)
    j = 1;
  end
  if(nargin<3)
    k = Params.knotDistrK;
  end
  doCentral = 0;
  if(j<-1)
    doCentral = 1;
  end
  j = abs(j);


  subtr = 0;
  if(doCentral) %compute central moment
    subtr = expect_k(pvec,1,k);
  end
  
  % subtract mean if j>1:
  k = k-subtr;
  dy = k.^j;
  dy = repmat(dy,Params.npp,1);
  if(isempty(pvec))  %return integrating vector:
    EV = dy;
  else
    EV = dot(dy,pvec);
  end
  
