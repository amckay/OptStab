function y0 = interp_vspline(vs,x0,doChangeOfVars)
% default for doChangeOfVars is true

  if nargin < 3
      doChangeOfVars = true;
  end
  
  if doChangeOfVars
      xx0 = log(x0 + 1);
  else
      xx0 = x0;
  end

  y0 = initsize(xx0,vs.x,vs.y);
  %y0 = zeros(size(x0));
  
  n = length(vs.y);
  if(any(x0>=vs.x(n)))
    max(x0)
    vs.x(n)
    error('too big CoH');
  end
  
  
  
  iPos = lookup(vs.x,xx0,1);
  
  y0 = vs.y(iPos) + max(xx0-vs.x(iPos),0) .* vs.Slope(iPos);

