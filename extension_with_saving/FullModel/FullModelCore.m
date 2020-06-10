if ~isfield(Params,'parstart')
    Params.parstart  = guessHouseholdPolicyRules();
end
% load('parameters/par0.mat','par');
% Params.parstart  = par;
%% COMPUTE STEADY STATE:




% get these by running Calibrate.m
wage = 0.83119;
beta = 0.99419;
AssetSupply = 0.21873;

% initial guesses 
% lambda0 = 0.719915510753214;
dividend0 = 0.142145022679927;
% mathcalH0 = 0.836013751445525;

% Look up prices
[ ~, R0, mathcalH0, ~, lambda0, ~, ~, ~, ~ ] = GetVarFromTable( beta, Params.b, false );


%%

Params.ststPrices = [lambda0; mathcalH0; dividend0];
% R0 = 1.005810416089240;
Params.fjacinv = diag([-1.0, -1.0, -1.0]);

R = fzero(@CheckR,[R0-0.0005, min(R0+0.0005,1/max(Params.beta)-0.0002)],[],wage,beta,AssetSupply)

optset('broyden','initb',Params.fjacinv)
optset('broyden','tol',1e-7);
[Params.ststPrices, ~, ststflag, ~, Params.fjacinv] = broyden(@CheckSteadyState, Params.ststPrices, 'counterfactual', R, wage, beta);
[ststres2,~,Env, X] = CheckSteadyState(Params.ststPrices, 'counterfactual', R, wage, beta);


assert(all(abs(ststres2) < 1e-4))
disp('Done with steady state')

wage = Env.xaggstst(strcmp(Params.aggnames,'wage'));
R = Env.xaggstst(strcmp(Params.aggnames,'R'));
lambda = Env.xaggstst(strcmp(Params.aggnames,'lambda'));
dividend = Env.xaggstst(strcmp(Params.aggnames,'dividend'));
B = Env.xaggstst(strcmp(Params.aggnames,'B'));
u = Env.xaggstst(strcmp(Params.aggnames,'u'));
mathcalH = Env.xaggstst(strcmp(Params.aggnames,'mathcalH')); 


%% Analyze the steady state

%open file for writing
fname = ['results/steady_state_' num2str(Params.b) '.txt'];
fid= fopen(fname,'w');


% expected marginal utility difference between unemployed and employed
[~, ~, ~, EMU_emp, EMU_unemp] = expect_C(Env.Dstst,Env.parstst,R,wage,dividend,lambda);
dispAndSave('expected marginal utility difference between employed and unemployed',fid)
dEMU = (EMU_unemp-EMU_emp)/EMU_emp;
dispAndSave(['Unweighted ' num2str(dEMU)],fid)

LiqAssetsToAnnualGDP = B/(4*(1-u)*mathcalH);
dispAndSave(['Liquid Assets To Annual GDP = ' num2str(LiqAssetsToAnnualGDP)],fid)
dispAndSave(['Steady state unemployment rate = ' num2str(u)],fid);

welfare =  expect_welfare(Env.Dstst,Env.parstst);
% welfare = welfare + Get_SteadyStateSkillWelfare( u ) / (1-Params.beta(end));
dispAndSave(['Steady state welfare = ' num2str(welfare)],fid);

% ***compute statistics of wealth distribution***
% compute median income
% m = get_medianIncome(D,par,wage,dividend,lambda,M);
% ShareLowWealth = expect_lowWealth(D,m);
% dispAndSave(['Fraction with less than one quarter''s median income = ' num2str(ShareLowWealth)],fid)

%household wealth is just bond wealth
[giniIndex, Lorenz] = gini(Env.Dstst, repmat(Params.knotDistrK,Params.npp,1), false);


dispAndSave(['Gini coefficient is ' num2str(giniIndex)], fid);
dispAndSave('Wealth quintiles:', fid)
dispAndSave('quintile    share of wealth', fid)
dispAndSave(['1           ' num2str(...
    Lorenz(find(Lorenz(:,1) > 0.2,1,'first'),2) - 0 ...
    )], fid)
dispAndSave(['2           ' num2str(...
    Lorenz((find(Lorenz(:,1) > 0.4,1,'first')),2) - Lorenz(find(Lorenz(:,1) <= 0.2,1,'last'),2)...
    )], fid)
dispAndSave(['3           ' num2str(...
    Lorenz(find(Lorenz(:,1) > 0.6,1,'first'),2) - Lorenz(find(Lorenz(:,1) <= 0.4,1,'last'),2)...
    )], fid)
dispAndSave(['4           ' num2str(...
    Lorenz(find(Lorenz(:,1) > 0.8,1,'first'),2) - Lorenz(find(Lorenz(:,1) <= 0.6,1,'last'),2)...
    )], fid)
dispAndSave(['5           ' num2str(...
    1 - Lorenz(find(Lorenz(:,1) <= 0.8,1,'last'),2)...
    )], fid)


% close steady state results file.
fclose(fid);
%% Linearization


Params.ntotal = length(X);

% the number of expectational errors
neta = (Params.nc+Params.nv)*Params.npp + Params.nAggEta;  

[X2,varindx] = setix(X,neta,Params.nz);

Env.varix = varindx;
[residStst,Env] = equ(X2,Env);
if any(abs(residStst) >1e-6)
    error('wrong stst');
end


njac = 20*ones(size(Env.varix.x));
njac(Env.iVarStatic)=500;
Env = dojacob(@equ,X2,Env,njac,Env);

%% Solution Linear RatEx Model 

disp('Solving the model exactly')
[G1,impact,eu,abseig] = solveexact(Env);
if(~all(eu==1))  %signal failure:
    G1 = [];
    error('sims:fails','sims failed');
end
    

%% Create observation matrix for all aggregate variables 

iA = Env.iVarAggr;  
niA = length(iA);

Hagg = sparse(1:niA,iA,ones(1,niA),niA,Params.ntotal)';
   


%% compute variances through simulation
disp('Compute variances through simulation...')

[ser, shocks] = simulateSystem(G1,impact,Hagg',modelOptions.sim_T, Params.Sigma,Params.seed);
ser = ser +  Env.xaggstst * ones(1,modelOptions.sim_T);

y = log(ser(strcmp(Params.aggnames,'Y'),:));
std_logy = std(y)
u = ser(strcmp(Params.aggnames,'u'),:);
std_u = std(u)

B = Env.xaggstst(strcmp(Params.aggnames,'B'));
R = Env.xaggstst(strcmp(Params.aggnames,'R'));
h = Env.xaggstst(strcmp(Params.aggnames,'mathcalH'));
lambda = Env.xaggstst(strcmp(Params.aggnames,'lambda'));
u = Env.xaggstst(strcmp(Params.aggnames,'u'));
wage = Env.xaggstst(strcmp(Params.aggnames,'wage'));


SetVarInTable( beta, Params.b , B, R, h, lambda, u, wage, dEMU, std_logy,std_u,welfare)
