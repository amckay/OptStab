%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION TO COMPUTE TRANSITION PROBABILITY MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = lineartrans(kgrid,k0)
  n = length(kgrid);
  k0 = max(min(k0,kgrid(n)),kgrid(1));
  iPos = lookup(kgrid,k0,0);
  iPos = reshape(iPos,1,length(k0));  % make sure it is column vector
  iPos = min(iPos,n-1);
  pHigh = (k0-kgrid(iPos))./(kgrid(iPos+1)-kgrid(iPos));
  pHigh = reshape(pHigh,1,length(k0));  % make sure it is column vector
  S.iFr = reshape(repmat(1:n,2,1),2*n,1);
  S.iTo = reshape([iPos;iPos+1],2*n,1);
  S.Val = reshape([1-pHigh;pHigh],2*n,1);