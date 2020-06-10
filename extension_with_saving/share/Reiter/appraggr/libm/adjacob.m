function J = adjacob(func,X,iDeriv,nmax,varargin)
  nx = length(X);
  if(isempty(iDeriv))
    iDeriv = 1:nx;
  end;
  n = length(iDeriv);
  Ind = sparse(nx,n);
  Ind(iDeriv,:) = speye(n);
  ndone = 0;
  J = [];
  % fprintf('adjacob: %d to do\n',n);
  while(ndone<n)
    if(isscalar(nmax))
      ii = ndone+1:min(n,ndone+nmax);
    else
      nmax0 = 1000;
      ii0 = ndone+1:min(n,ndone+nmax0);
      nCumul = cumsum(nmax(ii0));
      ii = ndone+(1:find(nCumul<=nmax0,1,'last'));
      assert(~isempty(ii));
    end
    xx = deriv1(X,[],Ind(:,ii));
    cput = cputime;
    ff = func(xx,varargin{:});
    if(isa(ff,'deriv1'))
      j = getjac(ff);
    else
      j = zeros(size(ff,1),length(ii));
    end
    J = [J j];
    ndone = max(ii);
    % fprintf('%d done in %0.2f sec\n',ndone,cputime-cput);
  end
