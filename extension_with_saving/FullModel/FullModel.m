clear all; close all; clc;
setpath()

global modelOptions;

% Some options about the solution algorthim
modelOptions.IRF_Length = 100;
modelOptions.sim_T = 10000;

run('parameters/initparams_search');

FullModelCore

disp('...done')

