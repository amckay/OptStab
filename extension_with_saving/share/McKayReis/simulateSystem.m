function [ser, shocks] = simulateSystem(A,B,H,T,Sigma,seed)

nshock = size(B,2);

if nargin == 6
    if verLessThan('matlab', '7.12')
        warning('calcvar:oldmatlab','Setting random seeds in this way is not supported in Matlab releases before R2011a.  Seed will not be set.'); 
    else
        rng('default')
        rng(seed)
    end
        
end

%simulate to calculate variances 
ser = zeros(size(H,1),T);
shocks = randn(T,nshock)';  %the RNG populates the matrix column by column so forming the matrix this way and then transposing it means that we get the same shocks even if we leave off some columns.
state = zeros(size(B,1),1);
for t = 1:T
    ser(:,t) = H*state;
    state = A*state + B*Sigma*shocks(:,t);
end
