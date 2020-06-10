function z = eulerres_stst(par,beta,R,wage,dividend,lambda,M)

  
  z = eulerres(par,par,beta,R,R,wage,wage,dividend,dividend,lambda,lambda,M,M);
  
end
  