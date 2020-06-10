% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function c = invmargutil(mu)
  global Params;
  assert(all(mu>0),'mu not positive');
  c = mu.^(-1/Params.gam);

