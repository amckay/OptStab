clear all; close all; clc;

setpath()
global modelOptions;

% Some options about the solution algorthim
modelOptions.IRF_Length = 100;
modelOptions.sim_T = 10000;

run('parameters/initparams_search');

blist = 0.84:-0.02:0.64;
nblist = numel(blist);

for iibb = 1:nblist
    Params.b = blist(iibb);
    
    
    
    FullModelCore;
    
end
disp('...done')

