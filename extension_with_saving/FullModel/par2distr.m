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
    D = reshape(D,Params.ndstst,Params.npp);
    
    par = reshape([-sum(par); par],Params.ndstst,Params.npp);
    
    D = D + par;
    D = D(:);


    assert(abs(sum(D) - 1)<1e-10,'assert failed in par2distr: distribution does not integrate to one.');

    
end