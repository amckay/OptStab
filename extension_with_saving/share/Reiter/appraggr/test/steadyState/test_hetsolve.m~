
%{
takes a guess of parameters, and uses fzero to find equilibrium capital stock
%}

function [par,Kequ,D] = test_hetsolve(Kequ,parstart)
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

 
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION THAT HAS TO BE ZERO IN STEADY STATE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
takes a capital stock, computes wage and R, solves steady state
consumption problem, computes saving transition matrix, finds the long
run distribution, aggregates capital and returns resid between actual
and initial guess.
%}
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
  z = test_eulerres(par,par,R,wage,Transf);

