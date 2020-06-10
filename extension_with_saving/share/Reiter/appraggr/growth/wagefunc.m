% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function wage = wagefunc(K,z);
  global Params;
  %adjust for_ differences from L=1
  KLratio = K / Params.ymean;
  wage = z*(1-Params.alpha)*Params.A*(KLratio).^Params.alpha;

