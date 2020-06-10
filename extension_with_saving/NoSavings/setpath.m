function setpath()
% Prepare the matlab path to solve model X where X can be:
% FullModel, RepAgent or HandToMouth





fs = filesep;
p_ = fileparts(pwd());



%check that we are in the correct dir
%assert(strcmp(p_(find(p_==fs,1,'last')+1:end),'autostab_repo'), 'It appears your working directory is not autostab_repo.')



%cleanup path from previous runs
restoredefaultpath;



addpath(p_);





% libraries:
addpath([p_ '/share/MirandaFackler']);
addpath([p_ '/share/Reiter/appraggr']);
addpath([p_ '/share/Reiter/appraggr/libaa']);
addpath([p_ '/share/Reiter/appraggr/rrqr']);
addpath([p_ '/share/Reiter/appraggr/libm']);
addpath([p_ '/share/Lengwiler']);
addpath([p_ '/share/dynare']);
addpath([p_ '/share/McKayReis']);
addpath([p_ '/NoSavings' ]);
clear aadir proj;


    
% to allow loading of gramian objects
a = gramian(0.5,0.5,100,1e-14); clear a;  



end