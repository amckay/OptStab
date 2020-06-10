function par = par2long(par)



global Params modelOptions;

nn = Params.nn;
nc = Params.nc;
nv = Params.nv;
npp = Params.npp;


[m n] = size(par);

if (m == (nn+nc+nv) * npp && n == 1)
    return
end

assert(m == nn+nc+nv && n == npp);

par = [reshape(par(1:nv,:),npp*nv,1);
    reshape(par(nv+1:nv+nn,:),npp*nn,1);
    reshape(par(nv+nn+1:end,:),npp*nc,1)];


end
