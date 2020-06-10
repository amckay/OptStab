
% Computes the labor supply conditional on being employed and computes the
% average search effort
% Inputs:
%   pvec:  vector of probabilities (histogram weights)
%   par: policy rule parameters

function H = expect_H(pvec,par)
global Params;

H = 0;
denom =0;

for ip = 1:Params.npp
    if Params.employed(ip)
        ind = (ip-1)*Params.ndstst+(1:Params.ndstst);


        N = nspline(par(Params.par_nind,Params.Nind(ip)));
        nsup = interp_nspline(N,Params.knotDistrK,true);

        H = H + sum(pvec(ind) .* nsup)*Params.Skill(ip);
        denom = denom + sum(pvec(ind));
    end
end
H = H / denom;



end
        
       

