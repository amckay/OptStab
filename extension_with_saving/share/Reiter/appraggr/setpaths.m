% Setting paths:

% set here the directory into which appraggr.tar.gz is installed:
aadir = '~/automatic_stabilization/autostab_repo/share/Reiter/appraggr/';

addpath(aadir);
% growth model:
addpath([aadir 'growth']);
% libraries:
addpath([aadir 'libaa']);
addpath([aadir 'rrqr']);
addpath([aadir 'libm']);
clear aadir;


% add directory where compecon toolbox (Miranda/Facler) is installed
%if all(computer() ~= 'MACI64')
%    addpath('/usr/local/matlabR2010a/toolbox/compecon/CEtools/');
%end
    
% to allow loading of gramian objects
a = gramian(0.5,0.5,100,1e-14); clear a;  

