global Params;


%% Scalars

%preferences

% Params.beta_dispersion = [0 0];
% Params.beta = 0.994 + Params.beta_dispersion;
% Params.beta = [0.994 0.994];
% Params.betaSwitch = 0;
% Params.nbeta = 1;
% Params.Qind = [2 2];
% Params.Nind = [1 1];

% 
% Params.beta_dispersion = [-0.0035 -0.0035 0 0];  %-0.02  yields Lorenz of 0.0083. Looking for .123
%                                                %-0.012 yields Lorenz of 0.0144
%                                                %-0.010 yields 0.0299
%                                                %-0.008 yields 0.0465
%                                                %-0.006        0.0649
%                                                %-0.004        0.1074
%                                                %-0.0035       0.1217
% Params.beta = 0.994 + Params.beta_dispersion;
% x   = 1/200*0.2/0.8;
% Params.betaSwitch = [1-x x;
%     1/200 1-1/200];
% clear x;
% Params.nbeta = 2;
% Params.Qind = [2 2 4 4];
% Params.Nind = [1 1 3 3];
% 
% 
% 




Params.beta_dispersion = [ -0.0002 -0.0002 -0.001 -0.001 0 0];
Params.beta = 0.994 + Params.beta_dispersion;
x   = 1/48;
Params.betaSwitch = [1-x x 0;
    0 1-x x;
    x 0 1-x];
Params.betaGroupSizes = invdistr(Params.betaSwitch);
Params.betaSwitchExpect = [1-x x 0;
    0 1-x x;
    0 0   1];
clear x;
Params.nbeta = 3;
Params.Qind = [2 2 4 4 6 6];
Params.Nind = [1 1 3 3 5 5];
Params.Skill = kron([0.808510632	1.063829737	1.063829737],[1 1]);







Params.sigma = 1;     %risk aversion
Params.gamma = 2; % labor supply parameter
Params.psy = 0.23573526400687425; % cost of unemployment
Params.kappa = 20.992453917558496;  % curvature of disutility of search effort


%markets
Params.mu = 1.2;
Params.theta = 0.2857142857142857;
Params.upsilon = 0.094148298362607782;
Params.ubar_groups = [0.077941213,0.04569041,0.059568377];
qM = 0.5914334174945006;
Params.upsilon_groups = kron(Params.ubar_groups .* qM ./(1-qM - Params.ubar_groups .* (1-qM)),[1 1]);
Params.psi_1 = 0.03088; 
Params.psi_2 = 1.0;
Params.zeta = 1.675977653631285;
Params.B = 0.076;  % assets to annual income

%government
OS_V1_Results;
Params.tau = taustar;
Params.b = bstar;
Params.chi = 0.2617792205155144;
Params.omegapi = 1.6578089766;
Params.omegau = 0.13321923864;


Params.rhoA = 0.9;
Params.rhoG = 0.9;
Params.rhoMP = 0.9;

% shock ordering A, MP, G
shock_covar = [  2.116e-05,   0.00000000e+00,   0.00000000e+00;
                  0.00000000e+00,   2.3965e-05,   0.00000000e+00;
                  0.00000000e+00,   0.00000000e+00,   2.14196937e-03];

Params.Sigma = chol(shock_covar);
Params.seed = 823672;
%% Household heterogeneity and income process

Params.npp = 2*Params.nbeta; %number of discrete household types.
                %you can use this to control "big" model with npp = 9
                %or small model (no skill differences) with npp = 3
Params.employed = repmat([true false],1,Params.nbeta);         

%% Approximation of policy rules and asset positions

Params.nc = 100; %number of points in approximation of individual cons func
Params.nn = 100; %number of points in approximation of individual labor supply func
Params.nv = 100; %number of points in approximation of individual labor supply func
Params.ncnv = Params.nc + Params.nn + Params.nv;

Params.ndstst = 250; %number of points in wealth histogram


Params.bhmin = 0;  %min bond position for household
Params.bhmax = 15;  %max bond position for household





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quadrature grid for_ wealth distribution:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.nStatesDistr = Params.ndstst*Params.npp;
% number knot points:
ndk = Params.ndstst ;
[Params.knotDistrK,Params.logshift] = makeknotd(Params.bhmin,Params.bhmax,ndk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KNOT POINTS FOR_ SAVINGS POLYNOMIAL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Params.xmin = 0.001;
Params.xmax = 0.9*Params.bhmax; % beyond that point, linear extrapolation is fine!

n1 = ceil(Params.nc/3);
x1 = linspace(Params.xmin,0.2,n1+1)';
x2 = logspaceshift(0.2,Params.xmax,Params.nc-n1,Params.logshift)';
Params.knotXi = [x1(1:end-1);x2];
%  Params.knotXi = logspaceshift(Params.xmin,Params.xmax,Params.nc,Params.logshift)';
Params.knotXi = Params.knotXi(2:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KNOT POINTS FOR_ LABOR SUPPLY POLYNOMIAL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Params.nxmin = 0.001;
Params.nxmax = Params.bhmax;
n1 = ceil(Params.nn/3);
x1 = linspace(Params.nxmin,0.2,n1+1)';
x2 = logspaceshift(0.2,Params.nxmax,Params.nn-n1,Params.logshift)';
Params.ngrid = [x1(1:end-1);x2];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KNOT POINTS FOR_ VALUE POLYNOMIAL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we define the value function as piece-wise linear in bhat = log(b+1)

Params.vxmin = log(1.001);
Params.vxmax = log(1+Params.bhmax); % beyond that point, log-linear extrapolation is fine!

Params.vgrid = linspace(Params.vxmin,Params.vxmax,Params.nv)';



%% Aggregate variables


Params.nz = 3; % number of exogenous var

Params.nAggEta = 2; % number of expectational errors in aggregate equations


exog = {'A', 'Gshock', 'mpshock'};
other_states = {};
static = {'C','u','I','ppi','Y','G','mathcalH','wage','pstar','M','J','B','dividend','lambda','R'};
dec = { 'pbarA', 'pbarB'};

Params.aggnames = {exog{:}, other_states{:}, static{:}, dec{:}};


Params.nagg = length(Params.aggnames);


% create indices
n = length(exog)+length(other_states);
Params.aggInd.state = logical([ones(1,n) zeros(1,Params.nagg-n)]);

n = length(dec);
Params.aggInd.dec = logical([zeros(1,Params.nagg-n) ones(1,n)]);

n = length(exog);
Params.aggInd.exog = logical([ones(1,n) zeros(1,Params.nagg-n)]);

Params.aggInd.static = ~logical(Params.aggInd.state + Params.aggInd.dec);



%% indexing

Params.par_sind = Params.nv+Params.nn+1:Params.nv+Params.nn+Params.nc;
Params.par_nind = Params.nv+1:Params.nv+Params.nn;
Params.par_vind = 1:Params.nv;