function [Cpath, Npath, Bpath, taxpath, U, NLF] = simulateforward(D0,parPath,ppi,wage,R,labormarketstatus)
%[Cpath, Npath] = simulateforward(parPath,ppi,wage,R,labormarketstatus)
%Simulates a distribution of households and computes total consumption,
%labor supply and asset position.
%
% Inputs: D0  -- initial distribution vector
% parPath -- output of solveback
% other -- see solveback for description
%
% Outputs: [X]path -- path for aggregate variable X over transition.


%given today's par and D, compute C, N, B
%find tomorrow's D

global Params;

T = size(parPath,2);

Cpath = zeros(1,T-1);
Npath = Cpath;
Bpath = Cpath;
taxpath = Cpath;
U = Cpath;
NLF = Cpath;

D = D0;

initialR = R(1);
R = R(2:end);

disp('simulating forward: ' )
for t = 1:T-1
    if t > 1
        InterestPaidOnAssets = R(t-1);
    else
        InterestPaidOnAssets = initialR;
    end
    [ Dprime, Cpath(t), Npath(t), Bpath(t), taxpath(t) ] = simulatestep( D,par2wide(parPath(:,t)),ppi(t),wage(t),InterestPaidOnAssets,labormarketstatus(t));
    D = Dprime;
    disc = sum(reshape(D,Params.ndstst,Params.npp));
    U(t) = sum(disc(Params.estatus==2));
    NLF(t) = sum(disc(Params.estatus==3));
end
    






