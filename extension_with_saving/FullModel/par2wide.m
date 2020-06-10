function par = par2wide(par)



global Params;

nn = Params.nn;
nc = Params.nc;
nv = Params.nv;
npp = Params.npp;


[m n] = size(par);


if m == nn+nc+nv && n == npp
    return
end




assert(m == (nn+nc+nv) * npp && n == 1);

par = [reshape(par(1:npp*nv),nv,npp);  % vpar in wide form
    reshape(par(npp*nv+1:npp*(nv+nn)),nn,npp);  % npar in wide form
    reshape(par(npp*(nv+nn)+1:end),nc,npp)];  % cpar in wide form
end
