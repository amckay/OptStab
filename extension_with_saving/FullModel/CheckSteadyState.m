function [res, B, Env, X] = CheckSteadyState(X,Mode,varargin)
% Checks market clearing for given prices
%
% Mode = counterfactual
% Takes a vector of R, lambda, H, Q and returns residuals of
% government budget
% average hours = H
% average search effort = Q
%
% varargin = wage, beta
%
%
% Mode = calibration
% Takes a vector of beta1, lambda, H, Q, wage and returns residuals of
% unemployment = utarget
% government budget clearing
% average hours = H
% average search effort = Q
%
% varargin = R, utarget
%

global Params;

% unpack

if strcmp(Mode , 'counterfactual')
    lambda = X(1);
    mathcalH = X(2);
    dividend = X(3);
    R = varargin{1};
    wage = varargin{2};
    beta_position = varargin{3};
    beta =  beta_position + Params.beta_dispersion;
elseif strcmp(Mode, 'calibration')
    lambda = X(1);
    mathcalH = X(2);
    dividend = X(3);
    M = X(4);
    R = varargin{1};
    utarget = varargin{2};
    beta_position = varargin{3};
    beta = beta_position + Params.beta_dispersion;

    wage = 1/Params.mu  - Params.upsilon * Params.psi_1 *M^Params.psi_2 / mathcalH;
end

disp(['lambda = ' num2str(lambda)])
disp(['mathcalH = ' num2str(mathcalH)])
disp(['dividend = ' num2str(dividend)])
disp(['R = ' num2str(R)])
disp(['wage = ' num2str(wage)])
disp(['beta = ' num2str(beta_position)])


if strcmp(Mode , 'counterfactual')
    M = ((1/Params.mu  - wage)*mathcalH/Params.psi_1/Params.upsilon)^(1/Params.psi_2);
end

[par,check] = broydn(@eulerres_stst,Params.parstart,[1e-11,0,1],beta,R,wage,dividend,lambda,M);
if(check~=0)
    %save broyError par;
    warning('compute_stst:broyerror','broydn not converged');
end


Params.parstart = par;
par = par2wide(par);

Pi = forwardmat(par,M);
D = invdistr(Pi);

Pi = forwardmatlabor(par,M);
Dtilde = Pi * D;

H = expect_H(Dtilde,par);  % average hours per employed worker


empbygrp = sum(reshape(Dtilde,Params.ndstst,Params.npp));
u = sum(empbygrp(~Params.employed));


A = 1;
M2 = ((A/Params.mu  - wage)*H/Params.psi_1/Params.upsilon)^(1/Params.psi_2);

Hires = Params.upsilon*(1-u);
Y = A * H * (1-u);
J = Params.psi_1 * M2.^(Params.psi_2) * Hires;
avgskillemployed = dot(empbygrp(Params.employed),Params.Skill(Params.employed));
dividend2 = ((Y-J) - wage * H *(1-u))/avgskillemployed;


[posttax_emp, posttax_unemp, pretax] = expect_tax(Dtilde,par,wage,dividend2,lambda);
[C, EInc] = expect_C(Dtilde,par,R,wage,dividend2,lambda);
B = expect_k(D)
% assert(abs(C - EInc - (R-1)*B) < 1e-4)
G = Params.chi * C;
lambda2 = (pretax-G - (R-1) *B)/(posttax_emp/lambda + posttax_unemp/lambda);

res = [lambda2-lambda; H-mathcalH; dividend2- dividend];
if strcmp(Mode, 'calibration')
      res = [res;u-utarget];
end

disp(res')

if(nargout>=3)
    % aggregate resource constraint (should clear by Walras)
    testResCon = Y-J-G-C
    assert(abs(testResCon) < 1e-3 * pretax)


    xaggstst= agg_stst(C, u, R, Y, G, mathcalH, wage, M, J, B, dividend, lambda);


    Env.Dstst = D;
    Env.parstst = par;
    Env.xaggstst = xaggstst;
    Env.beta = beta;



    X = [zeros(Params.nStatesDistr-1,1);xaggstst;par2long(par)];
end



end

