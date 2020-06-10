function c = get_c(x,savings,N,Rlag,wage,dividend,lambda,incomeIndex)
%compute consump, labor supply, and taxes
%inputs:
% x assets (vector)
% savings (vector)
% N  labor spline or vector of labor supplies
% R, wage prices


c = Rlag * x - savings + get_laborIncome(x,N,wage,dividend,lambda,incomeIndex);

