
% Computes the labor supply conditional on being employed and computes the
% average search effort
% Inputs:
%   pvec:  vector of probabilities (histogram weights)
%   par: policy rule parameters

function AggQ = expect_Q(pvec,par)
global Params;


AggQ = 0;
denom =0;
for ip = 1:Params.npp
        ind = (ip-1)*Params.ndstst+(1:Params.ndstst);
        
        if Params.employed(ip)
            probsearch = Params.upsilon;
        else
            probsearch = 1.0;
        end
        
        Q = nspline(par(Params.par_nind,Params.Qind(ip)));
        q = interp_nspline(Q,Params.knotDistrK,true);

        AggQ = AggQ + sum(pvec(ind) .* probsearch .* q);
        denom = denom + sum(pvec(ind) * probsearch);
    
end
AggQ = AggQ / denom;





end
        
       

