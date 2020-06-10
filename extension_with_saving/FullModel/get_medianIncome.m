function m = get_medianIncome(pvec,par,wage,dividend,lambda,M)


global Params;



npp = Params.npp;
par = par2wide(par);


Inc = zeros(Params.ndstst,npp);
Prob = zeros(Params.ndstst,npp);

%loop over discrete grid
for ib = 1:Params.nbeta
    ind = (ib-1)*Params.ndstst+(1:Params.ndstst);
    for iemp=1:2
    
        ip = (ib-1)*2+iemp;
    
        Q = nspline(par(Params.par_nind,Params.Qind(ip)));
        q = interp_nspline(Q,Params.knotDistrK,true);
        N = nspline(par(Params.par_nind,Params.Nind(ip)));

        Inc(:,ip) = get_laborIncome(Params.knotDistrK,N,wage,dividend,lambda,ip);
   
        if Params.employed(ip)
            probOfip = 1-Params.upsilon + Params.upsilon * q * M;
        else
            probOfip = Params.upsilon - Params.upsilon * q * M;
        end
        
        Prob(:,ip) = pvec(ind).*probOfip;
    
    end
    
end

Inc = Inc(:);
Prob = Prob(:);
Inc = Inc(Prob > 1e-7);
Prob = Prob(Prob > 1e-7);
[Inc, I] = sort(Inc);
Prob(:) = cumsum(Prob(I))/sum(Prob);

m = interp1(Prob,Inc,0.5);





end
