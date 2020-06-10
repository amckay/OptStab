function xout=subsref(x,s)
  nx = size(x);
  switch s.type
    case '()'
      % xout.v = x.v(s.subs{:});
      xout.v = subsref(x.v,s);
      Is = reshape(1:prod(size(x.v)),size(x.v));
      ii = subsref(Is,s);
      xout.d = x.d(ii(:),:);
    otherwise
      error('wrong index type for_ deriv1');
  end

  if(nnz(xout.d)==0)
    xout = xout.v;
  else
    xout = class(xout,'deriv1');
  end
  
