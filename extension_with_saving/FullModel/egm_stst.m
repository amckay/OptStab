function par = egm_stst(par,beta,R,wage,dividend,lambda,M,T)
  %compute steady state savings and labor supply using endog grid method
  
par = par2wide(par);
C1 = 0;
V1 = 0;

for t = 1:T
    [par, C2, V2] = egm(par,beta,R,wage,dividend,lambda,M,t > 100);
    testC = max(abs(C2(:)-C1(:))); 
    testV = max(abs(V2(:)-V1(:))); 
    if mod(t,100)==0
        disp(['Iteration ' num2str(t) ' -- ' num2str(testC) ' -- ' num2str(testV)])
    end
    C1 = C2;
    V1 = V2;
    if (testC < 1e-8) & (testV < 1e-8), break; end
end






    
end

