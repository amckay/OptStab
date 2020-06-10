function y = get_laborIncome(x,N,wage,dividend,lambda,incomeIndex)
%inputs:
% x  savings
% N  labor spline or vector of labor supplies
% R, wage prices

global Params;

if isstruct(N)
    n = interp_nspline(N,x,true); %true -> take max of zero
else
    n = max(N,0);
end

y =  (1- (1-Params.b) * ~Params.employed(incomeIndex)) * lambda *  (wage*n + dividend).^(1-Params.tau);

