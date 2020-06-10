% Given du/dc and n, this file computes c.
% Created Jan 20, 2012 by AGM
function c = invmargutilC(mu)
  global Params;
  assert(all(mu(:)>0),'mu not positive');
  
  
  
  c = mu.^(-1/Params.sigma);
  
