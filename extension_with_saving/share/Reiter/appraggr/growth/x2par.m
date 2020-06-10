function [D,par,dz,dtau,Kaggr,Env] = x2par(x,Env)
  global Params;
  nd = Params.nStatesDistr;
  nc = Params.nc*Params.npp;
  if(Params.iKS)
    D = x(1:nd);
  else
    % assert(Params.npp==1);
    D = par2distr(x(1:nd),Env);
  end


  Env.iExog = [nd+1 nd+2];
  dz = x(Env.iExog(1));
  dtau = x(Env.iExog(2));

  if(Params.iIncludeKAggr)
    Env.iVarStatic = [nd+3];
    Kaggr = x(Env.iVarStatic);
  else
    Env.iVarStatic = [];
    if(Params.iKS)
      % Kaggr = D(1)*Env.norm_h1;
      Kaggr = Env.mom2K(1:nd)'*x(1:nd);
    else
      Kaggr = make_h(1)'*x(1:nd) + Env.Kstst;
    end
  end
  Env.iVarAggr = [Env.iExog Env.iVarStatic];
  
  nx = length(x);
  Env.iVarDec = (nd+Params.nz+length(Env.iVarStatic)+1):nx;
  par = reshape(x(Env.iVarDec),Params.nc,Params.npp);
  Env.iVarState = [1:nd Env.iExog];
  Env.iVarBWS = [Env.iVarState Env.iVarStatic];  % backward and static

  %nb = length(Env.iVarBWS);
  %Params.HStateVars = eyeii(length(Env.iVarBWS),[nb-1 nb]);
  %Params.HDecVars = eye(length(Env.iVarDec));
  