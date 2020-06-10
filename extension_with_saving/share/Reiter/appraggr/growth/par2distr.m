% Matlab file for_ paper 
% "Approximate and Almost-Exact Aggregation in Dynamic Stochastic Heterogeneous-Agent Models"
% Michael Reiter, IHS, December 2009
% Last update: -
% Feel free to use, copy and modify at your own risk;
%   this program comes with NO WARRANTY WHATSOEVER
%
function D =  par2distr(par,Env)
  global Params;

  D = initsize(Env.Dstst,par);
  D(:) = Env.Dstst;
  D = reshape(D,Params.ndstst+1,Params.npp);
  par = reshape(par,Params.ndstst,Params.npp);
  D(1,:) = D(1,:)-sum(par);
  D(2:end,:) = D(2:end,:) + par;
  D = D(:);
  

