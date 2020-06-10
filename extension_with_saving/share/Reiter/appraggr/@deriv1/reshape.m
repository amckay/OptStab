function xOut = reshape(xIn,varargin)
  xOut.v = reshape(xIn.v,varargin{:});
  xOut.d = xIn.d;
  xOut=class(xOut,'deriv1');
