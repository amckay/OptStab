function y0 = interp_nspline(ns,x0,doMax)


  y0 = initsize(x0,ns.x,ns.y);
  %y0 = zeros(size(x0));
  
  n = length(ns.y);
  if(any(x0>=ns.x(n)))
    error('too big CoH');
  end
  
  
%   if (any(x0 < ns.x(1)))
%       disp(ns.x(1))
%       disp(min(x0))
%       error('Cannot extrapolate labor supply.');
%   end
  
  %iLow = x0 <= ns.x(1);
  %y0(iLow) = ns.y(1);
  
  %iPos = lookup(ns.x,x0(~iLow),0);
  iPos = lookup(ns.x,x0,1);
  %y0(~iLow) = ns.y(iPos) + (x0(~iLow)-ns.x(iPos)) .* ns.Slope(iPos);
  y0 = ns.y(iPos) + max(x0-ns.x(iPos),0) .* ns.Slope(iPos);
 
  if doMax
      y0 = max(y0,0);
  end