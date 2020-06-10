% Given c and n, this file computes du/dc.
% Created Jan 20, 2012 by AGM
function mu = margutilC(c,order)
  global Params;
  assert(all(c(:)>0),'c not positive');
  
  if ~exist('order','var') || order == 1
    mu = c.^(-Params.sigma);
  elseif order == 2
      mu = (-Params.sigma).*c.^(-Params.sigma-1);
  else
      error('unrecognized order passed to margutilC');
  end
      
  
