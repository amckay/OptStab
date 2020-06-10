%% Reiter accuracy check


% exact vs disaggregate (all aggregate variables)
%This may fail due to zero shock variance of G shocks
% as it is just a test, we can set that to something small instead of zero
Sigma = Params.Sigma;
Sigma(2,2) = 1e-4;

errAggr.diffIR_ExDis ...
    = evalmod(G1,impact,Hagg,{A_Dis},{B_Dis},{Hagg(:,1:niS)},[],Tmax,...
    Sigma,Params.betae);

Sigma = Params.Sigma;

%% Plot Exact and MPA solution for aggregate variables.  

myIRF_T = 500;

irEX = ir_sims(G1,impact,myIRF_T-1,Hagg,Params.SigmaEpsz); %exact
irEX = irEX{1};

irDis = ir_sims(A_Dis,B_Dis,myIRF_T-1,Hagg(:,1:niS),Params.SigmaEpsz); %MPA_disaggregated
irDis = irDis{1};

nplots = size(Hagg,1);

figure;
for i = 1:nplots
    subplot(ceil(nplots/4),4,i);
    y = irDis(:,i);
    y = [irEX(:,i) y]; %#ok
    plot(1:myIRF_T,y)
    title(Params.aggnames{i});
    if i == 1, legend('exact','model reduction'); end
end
