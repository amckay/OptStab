
% Computes the employment rate
% Inputs:
%   pvec:  vector of probabilities (histogram weights)

function emp = expect_employed(pvec)
global Params;

emp = 0;


for ip = 1:Params.npp
    if Params.employed(ip)
        ind = (ip-1)*Params.ndstst+(1:Params.ndstst);
        emp = emp + sum(pvec(ind));

    end
end
emp = emp / sum(pvec);



end
        
       

