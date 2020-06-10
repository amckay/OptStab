function [resid,Env] = equ(X,Env)

global Params;

%unpack
xcurr = X(Env.varix.x);
xlag = X(Env.varix.xlag);
eta = X(Env.varix.eta);
eps = X(Env.varix.eps);

[D,par,dz,xagg,Env] = x2par(xcurr,Env); %#ok
[Dlag,parlag,dzlag,xagglag] = x2par(xlag,Env); %#ok


% unpack stst, xagg and xagglag into local namespace
for nm_i = 1:length(Params.aggnames)
    eval([Params.aggnames{nm_i} '_stst = Env.xaggstst(nm_i);'])
    eval([Params.aggnames{nm_i} ' = xagg(nm_i);'])
    eval([Params.aggnames{nm_i} '_L = xagglag(nm_i);'])
end



nv = Params.nv*Params.npp;
nc = Params.nc*Params.npp;


%unpack expectational errors
eta_c_hhld = eta(1:nc);
eta_v_hhld = eta(nc+1:nc+nv);
eta_agg = eta(nc+nv+1:nc+nv+Params.nAggEta);


% auxiliary variables
%ET = expect_tax(D,par,wage,dividend,M);
% Q_impl = expect_Q(Dlag,par);

Pi = forwardmatlabor(par,M);
Dtilde = Pi * Dlag;
H_impl = expect_H(Dtilde,par);
u_impl = 1-expect_employed(Dtilde);

[posttax_emp, posttax_unemp, pretax] = expect_tax(Dtilde,par,wage,dividend,lambda);

empbygrp = sum(reshape(Dtilde,Params.ndstst,Params.npp));
avgskillemployed = sum(empbygrp(Params.employed).*Params.Skill(Params.employed));
% ---- aggregate equations (order them with backward, then static then forward) ---

resAgg = initsize(zeros(Params.nagg,1),X);

equCount = 0;
resAgg(equCount+1) = log(A)-Params.rhoA*log(A_L)-eps(1);  %AR(1)
resAgg(equCount+2) = log(mpshock) - Params.rhoMP*log(mpshock_L) -eps(2);
resAgg(equCount+3) = log(Gshock) - Params.rhoG*log(Gshock_L) -eps(3);
equCount = equCount + 3;


resAgg(equCount+1) = -wage + wage_stst*A*((1-u)/(1-u_stst)).^Params.zeta;
% resAgg(equCount+2 ) = -u + (u_L + Params.upsilon*(1-u_L))*(1-Q*M);
resAgg(equCount+2 ) = -u + u_impl;
resAgg(equCount+3 ) = -Y + A*mathcalH*(1-u);
resAgg(equCount+4 ) = -ppi + ((1-Params.theta)/(1-Params.theta*pstar.^(1/(1-Params.mu)))).^(1-Params.mu);  % pi
resAgg(equCount+5) = B_stst - B;  % check this
% resAgg(equCount+6) = C+G+J-Y;
resAgg(equCount+6) = C - expect_C(Dtilde,par,R,wage,dividend,lambda);
resAgg(equCount+7 ) = -G + Params.chi*C*Gshock;
resAgg(equCount+8 ) = -I + R_stst*ppi.^Params.omegapi*((1-u)/(1-u_stst)).^Params.omegau*mpshock;
resAgg(equCount+9 ) = -pstar + pbarA/pbarB; %pstar
resAgg(equCount+10) = -J + Params.psi_1 * M.^Params.psi_2 * (1-u - (1-Params.upsilon)*(1-u_L));
resAgg(equCount+11) = -pretax + G + (I_L/ppi-1)*B_stst + posttax_emp + posttax_unemp;  %government budget
resAgg(equCount+12) = -mathcalH + H_impl;
resAgg(equCount+13) = -dividend + ((Y-J) - wage*mathcalH*(1-u))/avgskillemployed;
resAgg(equCount+14) = B - expect_k(D);
resAgg(equCount+15) = R - I_L/ppi;
equCount = equCount + 15;

resAgg(equCount+1) = -pbarB_L + Y_L + (1-Params.theta)*ppi.^(-Params.mu/(1-Params.mu)) * ppi / I_L *pbarB +eta_agg(1); % pbarB    
resAgg(equCount+2) = -pbarA_L + Params.mu*Y_L*(wage_L * mathcalH_L+Params.psi_1*M_L.^Params.psi_2 - (1-Params.upsilon)*Params.psi_1*M.^Params.psi_2  )/ mathcalH_L / A_L + ...
    (1-Params.theta)*ppi.^(-Params.mu/(1-Params.mu)) * ppi/I_L*pbarA + eta_agg(2); % pbarA

equCount = equCount + 2;


% indexes
aggBW = 1:3;
aggStat = 4:19;
aggFW = 20:21;

Pi = forwardmat(par,M);
D2 = forward(Dlag,Pi);
resD = distr2par(D2,Env) - distr2par(D,Env);

%Rtmp = -Params.delta + z*Params.alpha*KL.^(-1 + Params.alpha);
%Rtmp_L = -Params.delta + z*Params.alpha*KL_L.^(-1 + Params.alpha);
%resCn = eulerres(parlag,par,1+Rtmp_L,1+Rtmp,wage_L,wage);

resVNC = eulerres(parlag,par,Env.beta,R_L,R,wage_L,wage,dividend_L,dividend,lambda_L,lambda,M_L,M);
    

% add expectational errors
tmp = Params.npp*Params.nn;
resVNC = resVNC + [eta_v_hhld; zeros(tmp,1);  eta_c_hhld];


resid = [resD;
    resAgg;
    resVNC];



% nbw = length(resD)+length(aggBW);
% nstat = length(aggStat);
% nfw = length(aggFW) + length(resVNC);
% 
% Env.iEquBW = 1:nbw;
% Env.iEquStat = nbw + 1 : nbw + nstat;
% Env.iEquFW = nbw + nstat + 1 : nbw + nstat + nfw;
% 
% Env.iEquBWS = [Env.iEquBW Env.iEquStat];
