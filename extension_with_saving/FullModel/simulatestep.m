function [ Dprime, C, N, Assets, taxpaid ] = simulatestep( D,par,ppi,wage,R,labormarketstatus)
%[ Dprime, C, N, Assets ] = simulatestep( D,par,ppi,wage,R,labormarketstatus)
%
% Inputs;
% D is distribution of assets entering date t before adjusting for inflation
% par is labor supply in t and savings from t to t+1
% ppi is inflation from t-1 to t
% wage is wage in t
% R is gross nominal interest from t-1 to t
% labormarketstatus describes transitions from t to t+1
%
% Outputs:
% Dprime is distribution entering t+1 before adjusting for inflation
% C is agg impatient cons in t
% N is agg impatient labor in t
% Assets is agg impatient saving from t to t+1
% tax payment is agg impatient income tax payment in t


global Params;


%rescale household savings by inflation
Dscaled = scaleassets(D,(1/ppi));




%household consumption
[C, ~, taxpaid] =  expect_C(Dscaled,par,R,wage);
C = Params.nu * C;
taxpaid = Params.nu * taxpaid;



N = expect_L(Dscaled,par,true)*Params.nu; % def of nh

Pi = forwardmat(1,par,labormarketstatus,true);
Dprime = forward(Dscaled,Pi);


Assets = expect_k(Dprime)*Params.nu;  % def of assetsh

assert(abs(sum(Dprime)-1) < 1e-6)

end

