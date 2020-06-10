function [c, n, margTax, taxpaid, taxableIncome] = get_cnt(x,savings,h,R,wage,incomeIndex,taxshock,dtau)
%compute consump, labor supply, and taxes
%inputs:
% x assets (vector)
% savings (vector)
% N  labor spline or vector of labor supplies
% R, wage prices

if ~exist('taxshock','var')
    taxshock = 0;
end

global Params;

if ~exist('dtau','var')
    dtau = 0;
end



if isstruct(h)
    n = interp_nspline(h,x,true); %true -> take max of zero
else
    n = max(h,0);
end

c = R * x - savings + (1- Params.b * ~Params.employed(incomeIndex)) * lambda *  (w*n + dividend)

effectiveWage = Params.skill(incomeIndex)*wage;

unemployed = ~all(Params.employed(incomeIndex));


taxableIncome = (R-1)*x + effectiveWage*n + Params.Tu(incomeIndex);
[margTax, taxpaid] = interp_tax(taxableIncome, Params.incometax,taxshock+dtau);
c = (x + taxableIncome - taxpaid - savings + Params.To(incomeIndex))/(1+Params.tauC+dtau);
