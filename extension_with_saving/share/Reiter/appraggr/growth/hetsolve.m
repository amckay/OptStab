% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function [par,Kequ,D,G1,impact,abseig,Env] = hetsolve(Kequ,parstart)
  global Params;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % FIND STEADY STATE:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(nargin>=2)  %Kequ and starting values for_ parameters given
    Params.parstart = parstart;
  else  %have to solve for_ Kequ
    %STARTING VALUES FOR_ CONSUMPTION FUNCTION:
    par0 = [0.1;0.95*Params.knotXi];
    Params.parstart = repmat(par0,Params.npp,1);
    if(nargin==1)
      Kstart = Kequ;
    else
      Kstart = [1;1.8];
    end
    % COMPUTE STEADY STATE:
    opts = optimset('fzero');
    opts.TolX = 1e-10;
    Kequ = fzero(@ststresid,Kstart,opts);
  end
  [resid,par,D,Env] = ststresid(Kequ);  % first use of Env!!;
  fprintf(1,'\nStSt done: Kequ = %0.6f\n',Kequ);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CHECK EULER RESIDUALS IN STEADY STATE:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % par(1) is kink point of consumption function_:
  xx2 = linspace(par(1),Params.xmax,10000)';
  R = 1+netintr(Kequ,1)-Params.taustst;
  wage = wagefunc(Kequ,1);
  Transf = Params.taustst*Kequ;
  resC = eulerres(par,par,R,wage,Transf,xx2);
  Params.eulerres = [max(abs(resC)) mean(abs(resC))];

  if(nargout<4)
    return;
  end

  % SOLVE LINEARIZED MODEL:


  if(Params.iIncludeKAggr)
    X = [zeros(Params.nStatesDistr,1);zeros(Params.nz,1);Kequ;par(:)];
  else
    X = [zeros(Params.nStatesDistr,1);zeros(Params.nz,1);par(:)];
  end
  ntotal = length(X);

  [X2,varindx] = setix(X,Params.nc*Params.npp,Params.nz); 
  
  Env.Kstst = Kequ;
  Env.varix = varindx;
  [residStst,Env] = equ(X2,Env);
  if(max(abs(residStst))>1e-6)
    error('wrong stst');
  end
  njac = 20*ones(size(Env.varix.x));
  njac(Env.iVarStatic)=500;
  Env = dojacob(@equ,X2,Env,njac,Env);
  G1 = []; impact = [];abseig = [];
  if(Params.npp==1)
    [G1,impact,eu,abseig] = solveexact(Env);
    if(~all(eu==1))  %signal failure:
      G1 = [];
      warning('sims failed');
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION THAT HAS TO BE ZERO IN STEADY STATE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resid,par,D,Env] = ststresid(K,Env)
  global Params;
  R = 1+netintr(K,1)-Params.taustst;
  wage = wagefunc(K,1);
  Transf = Params.taustst*K;

  % GIVEN K, SOLVE FOR THE CONSUMPTION FUNCTION BY COLLOCATION:
  % high precision in solution of savings function_ is required to achieve
  %   converge in outer loop!
  [par,check] = broydn(@eulerres_stst,Params.parstart,[1e-11,1,0],R,wage,Transf);
  if(check~=0)
    warning('broydn not converged');
  end
  par = reshape(par,Params.nc,Params.npp);
  Params.parstart = par;
  Pi = forwardmat(1,par,Transf,R,wage);
  D = invdistr(Pi);
  % D = invd2(Pi);
  clear Pi;

  K2 = expect_k(D);
  resid = K2 - K;
  fprintf(1,'K = %f, resid = %0.2e;  ',K-1, resid);

  if(nargout>3)
    Env.Dstst = D;
    Env.R = R;
    Env.wage = wage;
    Env.Transf = Transf;
    Env.parstst = par;
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION THAT DEFINES DISCRETE MODEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defines equations of stochastic model;
% Input:     X: all variables, current and lagged
% Outputs: residuals


function z = eulerres_stst(par,R,wage,Transf)
  z = eulerres(par,par,R,wage,Transf);

