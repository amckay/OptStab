function [res,Env, X] = CheckR(R,wage,beta,AssetSupply)
disp(sprintf('R %0.10f',R))
global Params;





optset('broyden','tol',1e-5);
optset('broyden','initb',Params.fjacinv);
[Params.ststPrices, ststres, ststflag, ~, Params.fjacinv] = broyden(@CheckSteadyState, Params.ststPrices, 'counterfactual', R, wage, beta);
[ststres2, B, Env, X] = CheckSteadyState(Params.ststPrices,'counterfactual', R, wage, beta);
assert(all(abs(ststres2) < 1e-4))



res = B-AssetSupply;
disp(sprintf('res %0.10f',res))
end

