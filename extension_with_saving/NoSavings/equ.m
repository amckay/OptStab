function [resid,Env] = equ(X,Env)

global Params;

%unpack
xcurr = X(Env.varix.x);
xlag = X(Env.varix.xlag);
eta = X(Env.varix.eta);
eps = X(Env.varix.eps);



% unpack stst, xagg and xagglag into local namespace
for nm_i = 1:length(Params.aggnames)
    eval([Params.aggnames{nm_i} '_stst = Env.aggstst(nm_i);'])
    eval([Params.aggnames{nm_i} ' = xcurr(nm_i);'])
    eval([Params.aggnames{nm_i} '_L = xlag(nm_i);'])
end





% ---- aggregate equations (order them with backward, then static then forward) ---

resAgg = initsize(zeros(Params.nagg,1),X);

equCount = 0;
resAgg(equCount+1) = log(A)-Params.rhoA*log(A_L)-eps(1);  %AR(1)
resAgg(equCount+2) = log(mpshock) - Params.rhoMP*log(mpshock_L) -eps(2);
resAgg(equCount+3) = log(Gshock) - Params.rhoG*log(Gshock_L) -eps(3);


equCount = equCount + 3;


resAgg(equCount+1 ) = -1/C_L +  Params.beta*(I_L/ppi/C  * (1 - u*(1-Params.b)) / (1-u_L*(1-Params.b)) * (1 - Params.upsilon*(1-Q*M)*(1-1/Params.b))  ) + eta(4); %Euler
resAgg(equCount+2 ) = -u + (u_L + (1-u_L)*Params.upsilon)*(1-Q*M);
resAgg(equCount+3 ) = -I + R_stst*ppi.^Params.omegapi*((1-u)/(1-u_stst)).^Params.omegau*mpshock;
resAgg(equCount+4 ) = -R + I / Eppi ;
resAgg(equCount+5 ) = -Rtrack1 + R_L;
resAgg(equCount+6 ) = -pstar + pbarA/pbarB; %pstar
resAgg(equCount+7 ) = -ppi + ((1-Params.theta)/(1-Params.theta*pstar.^(1/(1-Params.mu)))).^(1-Params.mu);  % pi
resAgg(equCount+8 ) = -Y + A*mathcalH*(1-u);
resAgg(equCount+9 ) = -G + Params.chi*C*Gshock;
resAgg(equCount+10) = -mathcalH + ((1-Params.tau)*wage* Y/A/(Y-J)).^(1/(1+Params.gamma));  % h dec rule
resAgg(equCount+11) = -Q.^Params.kappa + M*Vn;  % search dec rule
resAgg(equCount+12) = -wage + wage_stst*A*((1-u)/(1-u_stst)).^Params.zeta;
resAgg(equCount+13) = -H + 1-u - (1-Params.upsilon)*(1-u_L);
resAgg(equCount+14) = -J + Params.psi_1 * M.^Params.psi_2 * H;
resAgg(equCount+15) = C+G+J-Y;  %  agg res con
resAgg(equCount+16) = -dividend + (Y-J)/(1-u) - wage*mathcalH;
equCount = equCount + 16;

resAgg(equCount+1) = -pbarB_L + Y_L + (1-Params.theta)*ppi.^(-Params.mu/(1-Params.mu)) * ppi / I *pbarB +eta(1); % pbarB    
resAgg(equCount+2) = -pbarA_L + Params.mu*Y_L*((wage_L*mathcalH_L+Params.psi_1*M_L.^Params.psi_2-(1-Params.upsilon)*Params.psi_1*M.^Params.psi_2 )/ ( mathcalH_L * A_L) ) + ...
    (1-Params.theta)*ppi.^(-Params.mu/(1-Params.mu)) * ppi/I*pbarA + eta(2); % pbarA

resAgg(equCount+3) = Eppi_L - ppi + eta(3);
equCount = equCount + 3;


resAgg(equCount+1) = -Vn_L-log(Params.b) - mathcalH.^(1+Params.gamma)/(1+Params.gamma) + Params.psy + Params.beta*(1-Params.upsilon)*(1-Params.kappa/(1+Params.kappa)*Q*M)*Vn + eta(5);
equCount = equCount + 1;

resid =     resAgg;
    
