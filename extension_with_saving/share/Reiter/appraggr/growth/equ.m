function [resid,Env] = equ(X,Env)
  global Params;


  xcurr = X(Env.varix.x);
  xlag = X(Env.varix.xlag);
  eta = X(Env.varix.eta);
  eps = X(Env.varix.eps);

  [D,par,dz,dtau,Kaggr,Env] = x2par(xcurr,Env);
  [Dlag,parlag,dzlag,dtaulag,Kaggrlag] = x2par(xlag,Env);

  [R,wage,Transf] = factorprices(Kaggrlag,dz,dtau);

  resZ = [dz-Params.rhoz*dzlag-eps(1);
	  dtau-Params.rhot*dtaulag-eps(2)];

  if(Params.iKS)  %with Krusell/Smith method, dynamics of moments exogenously given:
    DKS = Env.Aks*([Dlag-Env.MomStst;dzlag;dtaulag]) + Env.Bks*eps;
    resD = -(D-Env.MomStst) + DKS(1:length(D));
  else
    Pi = forwardmat(0,par,Transf,R,wage);
    D2 = forward(Dlag,Pi);
    resD = distr2par(D2,Env) - distr2par(D,Env);
  end

  resC = eulerres(parlag,par,R,wage,Transf) + eta;
  if(Params.iIncludeKAggr)
    resK = Kaggr - expect_k(D);
  else
    resK = [];
  end
  resid = [resD;
	   resZ;
           resK;
	   resC];
    
  nbw = length(resD)+length(resZ)+length(resK);
  nc = length(resC);
  Env.iEquBW = 1:nbw;
  Env.iEquStat = [];
  Env.iEquFW = nbw+1:nbw+nc;

  Env.iEquBWS = [Env.iEquBW Env.iEquStat];
