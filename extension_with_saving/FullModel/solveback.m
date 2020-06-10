function parPath = solveback(parFinal,ppi,wage,R,labormarketstatus,betashock)
%parPath = solveback(parFinal,pricepath)
% SOLVES backwards from a period (e.g. steady state) in whcih there
% decision rules are given by parFinal and the prices are wage(end) and so
% on.  R(t) and ppi(t) are the gross nominal interest and inflation from
% date t-1 to t.  wage(t) is the wage in date t.  labormarketstatus(t) is the
% (log) labor market status from t-1 to t.  These four prices should be
% vectors with T elements.  parPath is then n by T where n is the length of
% the output of par2long(par).  wage(end) should be the wage in the steady
% state.

T = length(wage);

if ~exist('betashock','var')
    betashock = ones(size(R));
end


parPath = repmat(par2long(parFinal),1,T);

initialR = R(1);
R = R(2:end);

disp('solving back: ' )
for t = T-1:-1:1
    
    if mod(t,50) == 0
        disp(num2str(t))
    end
    
    if t < T-1
        ppinextnext = ppi(t+2);
    else
        ppinextnext = ppi(t+1);
    end
    
    if t > 1
        InterestPaidOnAssets = R(t-1);
    else
        InterestPaidOnAssets = initialR;
    end
    
    parPath(:,t) = par2long(egm(par2wide(parPath(:,t+1)),InterestPaidOnAssets,R(t),wage(t),wage(t+1),ppi(t+1),ppinextnext,labormarketstatus(t),betashock(t)));
    
    %FIXME: in Dynare, the timing of ii will probably differ from what I
    %am using here.
    
    
end
