%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION THAT COMPUTES MPC for a specific type of household
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   b:   current asset position
%   ip: current income type
%   trans: size of transfer to use for numerical calculation of marginal response
%   par:       parameter vector for_ this periods' savings function_
%   parnext:       parameter vector for_ next periods' savings function_
%   R,wage: 1+interest,wage paid this period.  i.e interest on assets entering the period
%   We assume steady state so R and wage are same next period
% outputs:
%   m: (1 x 2) mpc of dC/trans and then dC/dCashOnHand, where the latter
%       reflects changes in labor supply and taxes
%   check:  should be near zero, it is the residual in household budget
%           constraint
function [m, check] = mpc(b,ip, trans, par, parnext,R,wage)

global Params;

%use par to find initial conditios
par = par2wide(par);
S = savingspline(par(Params.par_sind,ip));
sthis = interp_savspline(S,b);

if Params.estatus(ip) == 1
    N = nspline(par(Params.par_nind,ip));
    [~, nthis] = get_cnt(b,sthis,N,R,wage,ip);
    X0 = [sthis;  nthis];
else
    X0 = sthis;
end


%compute c for (b, ip) unders standard policy rule (i.e. with trans = 0)
[X0,check] = broydn(@eulerres_mpc,X0,[1e-11,0,0],b,ip,0,parnext,R,wage);

%compute c for (b,ip) with trans extra after-tax income
[X1,check] = broydn(@eulerres_mpc,X0,[1e-11,0,0],b,ip,trans,parnext,R,wage);

%convert from savings to consumption
if Params.estatus(ip) == 1
    [c0, ~, tp0] = get_c_mpc(b,ip,0,X0(2), X0(1), R, wage);
    [c1, ~, tp1] = get_c_mpc(b,ip,trans,X1(2), X1(1), R, wage);
    
    effectiveWage = Params.skill(ip)*wage;
    denom = [ trans (trans + effectiveWage*(X1(2)-X0(2)) -(tp1-tp0))];
else
    [c0, ~, tp0] = get_c_mpc(b,ip,0,0, X0(1), R, wage);
    [c1, ~, tp1] = get_c_mpc(b,ip,trans,0, X1(1), R, wage);
    
    denom = [ trans (trans-(tp1-tp0))];
end


m = (c1-c0)./denom;
check = (X1(1) - X0(1) + (1+Params.tauC)*(c1-c0))/denom(2) - 1;

end

function res = eulerres_mpc(X,b,ip,trans, parnext,R,wage)

global Params modelOptions;

npp = Params.npp;

parnext = par2wide(parnext);


sthis = X(1);
if Params.estatus(ip) == 1
    nthis = X(2);
    res = initsize(zeros(2,1),X,b,ip,trans, parnext,R,wage);
else
    nthis = 0;
    res = initsize(zeros(1,1),X,b,ip,trans, parnext,R,wage);
end

[cthis, margTax] = get_c_mpc(b,ip,trans,nthis, sthis, R, wage);


if(any(cthis<0))  %signal inadmissible value to routine 'broydn';
    res = 1e100;
    return;
end




assets = sthis;

MUexp = 0;

laborMarketStatus = 0;

for jp=1:npp
    
    pp = transProb(ip,jp,laborMarketStatus, nthis);
    
    
    if(any(pp>0))
        Sn = savingspline(parnext(Params.par_sind,jp));
        Nn = nspline(parnext(Params.par_nind,jp));
        
        
        snext = interp_savspline(Sn,assets);
        
        
        [cnext, nnext, margtaxnext] = get_cnt(assets,snext,Nn,R,wage,jp);
        
        if(any(cnext<1e-8))  %signal inadmissible value to routine 'broydn';
            %           disp('B');
            res = 1e100;
            return;
        end
        
        if modelOptions.GHH
            MUnext = margutilC_GHH(cnext,nnext);
        else
            MUnext = margutilC(cnext);
        end
        
        MUexp = MUexp + pp.*MUnext .* (1+(1-margtaxnext)*(R-1));
        
        
    end
end

% Euler residual:
% expressed in relative consumption units:
%notice that R is incorporated into MUexp
if modelOptions.GHH
    res(1,1) = 1 - invmargutilC(Params.betah*MUexp)./(cthis - Params.psi1 .* nthis.^(1+Params.psi2)./(1+Params.psi2));
else
    res(1,1) = 1 - invmargutilC(Params.betah*MUexp)./cthis;
end


%calculations for labor supply residuals



if Params.estatus(ip) == 1
    
    if modelOptions.GHH
        laborResid = (1+Params.tauC) .* Params.psi1 .* nthis.^Params.psi2   ./(Params.skill(ip)*wage.*(1-margTax))-1;  %notice that n is allowed to be negative only here
    else
        laborResid = (1+Params.tauC) .* Params.psi1 .* nthis.^Params.psi2   ./(cthis.^(-Params.sigma).*Params.skill(ip)*wage.*(1-margTax))-1;  %notice that n is allowed to be negative only here
    end
    
    res(2,1) = laborResid;
    
end

end


function [cthis, margTax, taxpaid] = get_c_mpc(b,ip,trans,n, savings, R, wage)

global Params;

%use budget to find c
effectiveWage = Params.skill(ip)*wage;
nthis = max(n,0);
taxableIncome = (R-1)*b + effectiveWage*nthis + Params.Tu(ip);
[margTax, taxpaid] = interp_tax(taxableIncome, Params.incometax);
cthis = (b + taxableIncome - taxpaid - savings + Params.To(ip) ...
    + trans)/(1+Params.tauC);

end