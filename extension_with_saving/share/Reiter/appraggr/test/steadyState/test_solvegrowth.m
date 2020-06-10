
global Params;
npp = 1;  % permanent productivity states; npp=1 means only i.i.d shocks
           % npp=31 is the big calibration in the paper
sigtau = 0.01;  % stdev tax shock
freq = 4;  % frequency  (4 is quarterly etc.)
initparams(npp,freq,sigtau);

% COMPUTE STEADY STATE:
[par,Kequ,D] = hetsolve;

  

% COMPUTE DSGE
[par,Kequ,D,G1,impact,abseig,Env] = hetsolve(Kequ,par);
  

% index of static and backward looking variables:
iS = Env.iVarBWS;
% number of exogenous shocks:
nz = Params.nz;
n = length(iS);
% variable to be predicted: (aggregate capital)
Hpred = full(eyeii(n,Env.iVarStatic));  % aggr. capital is a
                                        % variable, with index Env.iVarStatic

