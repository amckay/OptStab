clear all; close all; clc;

global modelOptions;

% Some options about the solution algorthim
modelOptions.IRF_Length = 100;
modelOptions.sim_T = 50000;

run('parameters/initparams_search');
OS_V1_Results;  % utarget = 0.0641

Vn0 = 0.34;
h0 = 0.89;
q0 = 0.93;
lambda0 = 0.8;
M0 =  ((1/Params.mu  - Params.wbar)*h0/Params.psi_1/Params.upsilon)^(1/Params.psi_2);


ststPrices = [lambda0; h0; q0; Vn0; M0];    
fjacinv = diag([ -1.0, -1.0, -1.0, -1.0 -5.0]);

optset('broyden','initb',fjacinv)
[X, ststres, ststflag, ~, Params.fjacinv] = broyden(@compute_stst_v2, ststPrices, 'calibration', utarget);

lambda = X(1);
h = X(2); % average hours among employed
q = X(3); % average search effort among unemployed
Vn = X(4);
M = X(5);


wage = 1/Params.mu  - Params.upsilon * Params.psi_1 *M^Params.psi_2 / h;
disp(['wage = ' num2str(wage,12)])

