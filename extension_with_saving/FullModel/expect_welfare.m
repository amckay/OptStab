function welfare = expect_welfare(pvec,par)

global Params;


par = par2wide(par);


welfare = 0;

%loop over discrete grid
for ip = 1:Params.npp
    ind = (ip-1)*Params.ndstst+(1:Params.ndstst);
    
    V = vspline(par(Params.par_vind,ip));
    v = interp_vspline(V,Params.knotDistrK);

    welfare = welfare + sum(pvec(ind) .*  v);

    
end


end
