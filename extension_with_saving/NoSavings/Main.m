clear all; close all; clc;

global modelOptions;

setpath()

% Some options about the solution algorthim
modelOptions.IRF_Length = 100;
modelOptions.sim_T = 50000;

run('parameters/initparams_search');


blist = [0.7:0.01:0.9 0.783];
iibb = 10;
for iibb = 1:numel(blist)
    Params.b = blist(iibb);
%% COMPUTE STEADY STATE:


Vn0 = 0.34;
h0 = 0.89;
q0 = 0.93;
lambda0 = 0.8;

ststPrices = [lambda0; h0; q0; Vn0];    
fjacinv = diag([ -1.0, -1.0, -1.0, -1.0]);
optset('broyden','initb',fjacinv)
[ststPrices, ststres, ststflag, ~, Params.fjacinv] = broyden(@compute_stst_v2, ststPrices, 'counterfactual');

[ststres,xaggstst, Env, X, lambda, welfare] = compute_stst_v2(ststPrices,'counterfactual');

disp('Done with steady state')

u = xaggstst(strcmp(Params.aggnames,'u'));
welfare = welfare + Get_SteadyStateSkillWelfare( u ) / (1-Params.beta);
disp(['Steady state welfare = ' num2str(welfare)])

%% Linearization


Params.ntotal = length(X);

% the number of expectational errors
neta =  Params.nAggEta;  

[X2,varindx] = setix(X,neta,Params.nz);

Env.varix = varindx;
[residStst,Env] = equ(X2,Env);
if any(abs(residStst) >1e-6)
    error('wrong stst');
end


njac = 20*ones(size(Env.varix.x));
Env = dojacob(@equ,X2,Env,njac,Env);

%% Solution Linear RatEx Model -- first exactly

disp('Solving the model exactly')
[G1,impact,eu,abseig] = solveexact(Env);
if(~all(eu==1))  %signal failure:
    G1 = [];
    error('sims:fails','sims failed');
end
    




%% Create observation matrix for all aggregate variables 

Hagg = eye(Params.ntotal);
   


%% compute variances through simulation
disp('Compute variances through simulation...')

[ser, shocks] = simulateSystem(G1,impact,Hagg',modelOptions.sim_T, Params.Sigma,Params.seed);
ser = ser +  Env.aggstst * ones(1,modelOptions.sim_T);

y = log(ser(strcmp(Params.aggnames,'Y'),:));
u = ser(strcmp(Params.aggnames,'u'),:);
std_logy = std(y)
stdu = std(u)

B = 0;
R = xaggstst(strcmp(Params.aggnames,'R'));
h = xaggstst(strcmp(Params.aggnames,'mathcalH'));
q = xaggstst(strcmp(Params.aggnames,'Q'));
u = xaggstst(strcmp(Params.aggnames,'u'));
wage = xaggstst(strcmp(Params.aggnames,'wage'));
EMU = 1/Params.b -1;


welfareWithVol = welfare + Params.CostOfVolatility * std_logy

SetVarInTable( Params.beta(end), Params.b , B, R, h,  lambda, u, wage, EMU, std_logy,stdu,welfare)
 end
disp('...done')

