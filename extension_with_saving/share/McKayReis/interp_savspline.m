function y0 = interp_savspline(S,x0)
  y0 = initsize(x0,S.x,S.y);
  y0(x0<S.x(1)) = 0;  % for lower CoH, no saving
  n = length(S.y);
  if(any(x0>=S.x(n)))
    error('too big CoH');
  end
  iMid = x0>S.x(1);
  iPos = lookup(S.x,x0(iMid),0);
  y0(iMid) = S.y(iPos) + (x0(iMid)-S.x(iPos)) .* S.Slope(iPos);
