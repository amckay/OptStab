function [res] = CheckBeta(beta,R0,utarget)

global Params;


optset('broyden','initb',Params.fjacinv)
[Params.ststPrices, ststres, ststflag, ~, Params.fjacinv] = broyden(@CheckSteadyState, Params.ststPrices, 'calibration', R0, utarget, beta);


mathcalH = Params.ststPrices(2);
M = Params.ststPrices(4);


wage = 1/Params.mu  - Params.upsilon * Params.psi_1 *M^Params.psi_2 / mathcalH;


%--
par = par2wide(Params.parstart);
Pi = forwardmat(par,M);
D = invdistr(Pi);
D = reshape(D,Params.ndstst,Params.npp);

summer = zeros(6,3); summer(1:2,1) = 1; summer(3:4,2) = 1; summer(5:6,3) = 1;
D2 = D * summer;
D2  = bsxfun(@rdivide,cumsum(D2),sum(D2));

[~,i] = min(abs(D2-0.5));
medianAssets = Params.knotDistrK(i) / ( wage * (1-utarget) * mathcalH * 4);
res = medianAssets(3) - 0.123924275

end

