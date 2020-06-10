function [R,wage,Transf] = factorprices(Dlag,dz,dtau)
  global Params;

  if(length(Dlag)<Params.ndstst-2)  %with Krusell/Smith method, aggregate K is first moment:
    Klag = Dlag(1);
  else
    Klag = expect_k(Dlag);
  end
  assert(Klag>0);

  tau = Params.taustst+dtau;
  
  R = 1+netintr(Klag,1+dz)-tau; 
  wage = wagefunc(Klag,1+dz);
  Transf = tau*Klag;
  if(Params.iWageTax)
    wage = wage + Transf;
    Transf = 0;
  end
