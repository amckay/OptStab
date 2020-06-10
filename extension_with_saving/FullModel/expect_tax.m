
% Computes integrals in the government budget constraint
% Inputs:
%   pvec:  vector of probabilities (histogram weights)
%   par: policy rule parameters

function [posttax_emp, posttax_unemp, pretax] = expect_tax(pvec,par,wage,dividend,lambda)
global Params;


posttax_emp = 0;
posttax_unemp = 0;
pretax = 0;


for ip = 1:Params.npp
    ind = (ip-1)*Params.ndstst+(1:Params.ndstst);
    
    
    N = nspline(par(Params.par_nind,Params.Nind(ip)));
    nsup = interp_nspline(N,Params.knotDistrK,true);
        
    if Params.employed(ip)
        pretax = pretax + sum(pvec(ind) .* (wage*nsup+dividend)*Params.Skill(ip));
        posttax_emp = posttax_emp + sum(pvec(ind) .* lambda .* ((wage*nsup+dividend)*Params.Skill(ip)).^(1-Params.tau));
    else
        posttax_unemp = posttax_unemp + sum(pvec(ind) .* Params.b .* lambda .* ((wage*nsup+dividend)*Params.Skill(ip)).^(1-Params.tau));
    end
end

end