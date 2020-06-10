% creates inputs for file dosims.m:
% Outputs:
%   X2: vector with all variables (including shocks)
%   indx: index structure (where in the vector the variables are)
function [X2,indx] = setix(X,nEta,nEps)

  nx = length(X);
  indx.x = 1:nx;
  indx.xlag = nx+indx.x;
  nn = max(indx.xlag);
  indx.eta = nn+1:nn+nEta;
  indx.eps = nn+nEta+1:nn+nEta+nEps;

  nX2 = max([max(indx.x),max(indx.xlag),max(indx.eps),max(indx.eta)]);
  X2 = zeros(nX2,1);
  X2([indx.xlag indx.x]) = repmat(X,2,1);

