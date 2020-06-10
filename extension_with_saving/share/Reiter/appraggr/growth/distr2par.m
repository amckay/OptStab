% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function par =  distr2par(D,Env)
  global Params;

  Ddev = D - Env.Dstst;
  Ddev = reshape(Ddev,Params.ndstst+1,Params.npp);
  par = Ddev(2:end,:);
  par = par(:);
  
