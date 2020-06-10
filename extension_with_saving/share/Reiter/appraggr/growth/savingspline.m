% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
% Input: 
%   par: finite parameterization of the consumption function
% Outputs:  
%   x:   knot points of spline
%   s:   saving at knot points
%   s2:  second derivative of saving at knot points
% The ouputs define a spline approximation for the consumption (or savings) function
function S = savingspline(par,X)
  global Params;

  if(nargin==2)
    par = [X(1);par(2:end)];
    knotXi = X(2:end)-X(1);
  else
    knotXi = Params.knotXi;
  end

  S.xcrit = par(1);
  S.x = [S.xcrit ; S.xcrit+knotXi; 1e8];  %add end point 1e8/
  S.y = [0;par(2:end);0];  
  n = length(S.y);
  S.y(n) = S.y(n-1) + (S.x(n)-S.x(n-1))*(S.y(n-1)-S.y(n-2))/(S.x(n-1)-S.x(n-2));

  % S.iSlope = diff(S.x) ./ diff(S.y);
  S.Slope = diff(S.y) ./ diff(S.x);

