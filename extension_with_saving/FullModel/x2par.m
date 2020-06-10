function [D,par,dz,xaggr,Env] = x2par(x,Env)
  global Params;
  
  nd = Params.nStatesDistr-1;
  
  nc = Params.nc*Params.npp;
  nv = Params.nv*Params.npp;
  nn = Params.nn*Params.npp;
  
  
  D = par2distr(x(1:nd),Env);
  

  tmp_n = nd+Params.nagg;

  Env.iExog = nd + find(Params.aggInd.exog);
  Env.iVarStatic = nd + find(Params.aggInd.static);
  Env.iVarAggr = [nd+1:tmp_n];
  Env.iVarAggrBWS = [nd+1:tmp_n-length(find(Params.aggInd.dec))];
  
  xaggr = x(Env.iVarAggr);
  dz = x(Env.iExog);
 
  
  
  Env.iVarDec = [nd+find(Params.aggInd.dec) tmp_n+1:tmp_n+nn+nc+nv];
 
    %tmp_n = nd + Params.nagg;
  par = par2wide(x(tmp_n+1:end));
    
  
  Env.iVarState = [1:nd nd+find(Params.aggInd.state)];
  Env.iVarBWS = [Env.iVarState Env.iVarStatic];  % backward and static

  %nb = length(Env.iVarBWS);
  %Params.HStateVars = eyeii(length(Env.iVarBWS),[nb-1 nb]);
  %Params.HDecVars = eye(length(Env.iVarDec));
  