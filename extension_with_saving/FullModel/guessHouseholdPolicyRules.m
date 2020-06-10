function [ parstart] = guessHouseholdPolicyRules()


global Params;

beta = 0.994 + Params.beta_dispersion;
wage = 0.82;
R = 1.005;
dividend = 0.167;
M = 0.6537;
lambda = 0.7455;

%initial guess
pars = [0; Params.knotXi]; %first element is sup {a :a'(a) = 0}
parn = 0.31*ones(Params.nn,1);
parv = ones(Params.nv,1);

par0 = [parv; parn; pars];
par0 = repmat(par0,1,Params.npp);



if size(par0,2) == 1
    par0 = par2wide(par0);
end

%to get a good initial guess, use endog grid method but it does not need to
%converge fully
disp('Using endog. grid method to find a good initial guess for policy rules.');
parstart = egm_stst(par0,beta,R,wage,dividend,lambda,M,1500);

parstart = par2long(parstart);


end

