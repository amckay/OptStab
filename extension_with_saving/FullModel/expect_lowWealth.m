function EV = expect_lowWealth(pvec,thresh)
  global Params;
  
  dy = Params.knotDistrK;

  
  dy = repmat(dy,Params.nbeta,1);
  
  EV = sum(pvec.*(dy <= thresh));

  
