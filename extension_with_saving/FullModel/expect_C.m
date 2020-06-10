% [EC Emargtax Etaxpaid] = expect_C(pvec,par,R,wage)
% Computes the aggregate consumption and taxes of households
% Inputs:
%   pvec:  vector of probabilities (histogram weights)
%   par: policy rule parameters
%   R: interest rate
%   wage
%   disc: vector of discrete states levels over which to integrate
%
% Alisdair McKay 2/24/12

function [EC, EInc, ES, EMU_emp, EMU_unemp] = expect_C(pvec,par,R,wage,dividend,lambda)

global Params;



par = par2wide(par);


EC = 0;
EInc = 0;
ES = 0;
EMU_emp = 0;
EMU_unemp = 0;
denom_emp = 0;
denom_unemp = 0;

%loop over discrete grid
for ip = 1:Params.npp
    ind = (ip-1)*Params.ndstst+(1:Params.ndstst);
    
   
    N = nspline(par(Params.par_nind,Params.Nind(ip)));

    S = savingspline(par(Params.par_sind,ip));
    
    
    %we have end of period assets in last period and this period,
    %and we have to figure out this period's consumption:
    xthis = Params.knotDistrK;
    sthis = interp_savspline(S,xthis);
    c = get_c(xthis,sthis,N,R,wage,dividend,lambda,ip);

    li = get_laborIncome(xthis,N,wage,dividend,lambda,ip);
    
    EC = EC + sum(pvec(ind) .* c);
    if nargout >= 2
        EInc = EInc + sum(pvec(ind)  .*li);
        ES = ES + sum(pvec(ind)  .* sthis);
    end    
    if nargout >= 4
        if Params.employed(ip)
            EMU_emp = EMU_emp + sum(pvec(ind) .* 1./c);
            denom_emp = denom_emp + sum(pvec(ind));
        else
            EMU_unemp = EMU_unemp + sum(pvec(ind) .* 1./c);
            denom_unemp = denom_unemp + sum(pvec(ind));
        end
    end
end
if nargout >= 4
    EMU_emp = EMU_emp / denom_emp;
    EMU_unemp = EMU_unemp / denom_unemp;
end


