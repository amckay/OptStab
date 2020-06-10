global Params;


%% Scalars

%preferences
Params.beta = 0.9818;  

Params.sigma = 1;     %risk aversion
Params.gamma = 2; % labor supply parameter
Params.psy = 0.23573526400687425; % cost of unemployment
Params.kappa = 20.992453917558496;  % curvature of disutility of search effort


%markets
Params.mu = 1.2;
Params.theta = 0.2857142857142857;
Params.upsilon = 0.094148298362607782;
Params.psi_1 = 0.03088; 
Params.psi_2 = 1.0;
Params.zeta = 1.675977653631285;
Params.wbar = 0.831258437422;
                
%government
Params.debt2GDP = 1.69781919145;  % stst debt to GDP ratio
OS_V1_Results;
Params.b = bstar;
Params.tau = taustar;
Params.chi = 0.2617792205155144;
Params.Bbar = 0;
Params.omegapi = 1.6578089766;
Params.omegau = 0.13321923864;



Params.rhoA = 0.9;
Params.rhoG = 0.9;
Params.rhoMP = 0.9;

% shock ordering A, MP, G
shock_covar = [  2.116e-05,   0.00000000e+00,   0.00000000e+00;
                  0.00000000e+00,   1.9165e-05,   0.00000000e+00;
                  0.00000000e+00,   0.00000000e+00,   2.14196937e-03];

Params.Sigma = chol(shock_covar);
Params.seed = 823672;

Params.CostOfVolatility = 10;


%% Aggregate variables


Params.nz = 3; % number of exogenous var

Params.nAggEta = 5; % number of expectational errors in aggregate equations


exog = {'A', 'Gshock', 'mpshock'};
other_states = {};
static = {'u','I','R','Rtrack1','ppi','Y','G','mathcalH','wage','pstar','M','Q','H','J','dividend'};
dec = { 'C','pbarA', 'pbarB', 'Eppi', 'Vn'};

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

