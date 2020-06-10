% INITIALIZE PARAMETER VALUES
% sets global structure Params;
function initparams(npp,freq,sigtau);
global Params;

Params.npp = npp; 
Params.freq = freq; 
Params.SigmaEps = [0.007 0; 0 sigtau].^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tax in steady state:
Params.taustst = 0.0;
% borrowing constraint:
Params.kmin = 0;


Params.iWageTax = 1;
% Model parameters:
Params.alpha = 1/3;
Params.beta = 0.99^(4/Params.freq);
Params.delta = 0.1/Params.freq;
Params.gam = 1;
% number of shocks:
Params.nz = 2;
% correlation shocks:
Params.rhoz = 0.95^(4/Params.freq);
Params.rhot = 0.;
Params.rhozz = diag([Params.rhoz Params.rhot]);
Params.iKS = 0;  % will be nonzero for Krusell/Smith method
Params.iIncludeKAggr = 1;  % use Kaggr as extra variable

% income process:
caliinc(Params.npp);


% number of exogenous shocks:
Params.nz = 2;
Params.nc = 100;
Params.ndstst = 1000;


% GRIDS:

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % DISCRETE APPROXIMATION OF INCOME PROCESS:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % lognormal distribution:
  Params.ymean = 1;
  Params.ymax = max(Params.xY);
  Params.ymin = min(Params.xY);
  
  % Upper limit capital:
  if(Params.npp==1)
    Params.kmax = 10*Params.ymax;
  else
    Params.kmax = 100;  
  end

  % CHECK COMPATIBILITY OF R AND K
  % take into account that Y does not add up exactly to 1 because of truncation
  R = 1/Params.beta;
  Params.A = (R-1+Params.delta)/Params.alpha;
  KLratio = Params.ymean;
  assert(abs(R-1-netintr(KLratio,1))<1e-8);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % quadrature grid for_ wealth distribution:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Params.nStatesDistr = Params.ndstst*Params.npp;  %makes only sense for npp==1; otherwise_:DSF!
  % number knot points:
  ndk = Params.ndstst + 1;
  [Params.knotDistrK,Params.logshift] = makeknotd(Params.kmin,Params.kmax,ndk);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % KNOT POINTS FOR_ SAVINGS POLYNOMIAL:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Params.xmin = 0.001;
  Params.xmax = Params.kmax*0.5; % beyond that point, linear extrapolation is fine!
  n1 = ceil(Params.nc/3);
  x1 = linspace(Params.xmin,0.2,n1+1)';
  x2 = logspaceshift(0.2,Params.xmax,Params.nc-n1,Params.logshift)';
  Params.knotXi = [x1(1:end-1);x2];
  %  Params.knotXi = logspaceshift(Params.xmin,Params.xmax,Params.nc,Params.logshift)';
  Params.knotXi = Params.knotXi(2:end);

  % use polynomials in logs, for moments:
  Params.momtype = 'pl';