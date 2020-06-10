clear all; close all; clc;

run('parameters/initparams_search');

Params.parstart  = guessHouseholdPolicyRules();

%% COMPUTE STEADY STATE:
% OS_V1_Results; %utarget = 0.0641

utarget = 0.0652;  % idea is to look at the plot of underinsurance and try to get the
                   % savings case to roughly line up with the others.

R = 1.005;
beta = 0.9946;   

% lambda0 = 0.72;
% mathcalH0 = 0.85;
% Q0 = 0.92;
% R0 = 1.005;
% M0 = 0.6537;
% 
% 
% 
% Jinv = diag([ -2/3, -1.0, -1.0 -5.0]);
% optset('broyden','initb',Jinv)
% [ststPrices, ststres, ststflag, ~, fjacinv] = broyden(@CheckSteadyState, [lambda0; mathcalH0; Q0; M0], 'calibration', R0, utarget, beta);

 
Params.fjacinv = diag([ -2/3, -1.0, -1.0 -5.0]);
 
lambda0 = 0.72;
mathcalH0 = 0.85;
dividend0 = 0.12;
M0 = 0.6537;
Params.ststPrices = [lambda0; mathcalH0; dividend0; M0];


beta_position = fzero(@CheckBeta,0.99421+[-0.0005, 0.0003],[],R,utarget)
% beta_position = 0.99419;
beta = beta_position + Params.beta_dispersion;

optset('broyden','initb',Params.fjacinv)
[Params.ststPrices, ~, ststflag, ~, Params.fjacinv] = broyden(@CheckSteadyState, Params.ststPrices, 'calibration', R, utarget, beta_position);
[ststres,AssetSupply] = CheckSteadyState(Params.ststPrices, 'calibration', R, utarget, beta_position);
disp(['AssetSupply = ' num2str(AssetSupply)])

ststPrices = Params.ststPrices;

lambda = ststPrices(1);
mathcalH = ststPrices(2);
dividend = ststPrices(3);
M = ststPrices(4);

wage = 1/Params.mu  - Params.upsilon * Params.psi_1 *M^Params.psi_2 / mathcalH;
% disp(['beta = ' num2str(ststPrices(1))])
disp(['wage = ' num2str(wage)])


par = broydn(@eulerres_stst,Params.parstart,[1e-11,0,1],beta,R,wage,dividend,lambda,M);
par = par2wide(par);
Pi = forwardmat(par,M);
D = invdistr(Pi);
D = reshape(D,Params.ndstst,Params.npp);

summer = zeros(6,3); summer(1:2,1) = 1; summer(3:4,2) = 1; summer(5:6,3) = 1;
D2 = D * summer;
D2  = bsxfun(@rdivide,cumsum(D2),sum(D2));

[~,i] = min(abs(D2-0.5));
medianAssets = Params.knotDistrK(i) / ( wage * (1-utarget) * mathcalH * 4)


% expected marginal utility difference between unemployed and employed
[~, ~, ~, EMU_emp, EMU_unemp] = expect_C(D(:),Params.parstart,R,wage,dividend,lambda);
disp('expected marginal utility difference between employed and unemployed')
dEMU = (EMU_unemp-EMU_emp)/EMU_emp;
disp(['Unweighted ' num2str(dEMU)])
