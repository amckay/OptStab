function [xaggstst, welfare]= agg_stst(R,lambda,h,q,Vn,wage)

    global Params;
    
    [u, M, Y, dividend,J] = compute_steady_state_supply_side(h,q,wage);
    
    
    A = 1;
    Gshock = 1;
    mpshock = 1;
    
    C = lambda*(wage*h + dividend)^(1-Params.tau)*(1-u+u*Params.b);
    G = Params.chi * C;
    
    I = R;
    Rtrack1 = R;
    ppi = 1;
    Eppi = 1;
    mathcalH = h;
    pstar = 1;
    Q = q;
    H = Params.upsilon*(1-u); 
    
    pbarB = Y/(1-  (1-Params.theta)/I);  
    pbarA = pbarB;
    
    xaggstst = zeros(Params.nagg,1);
    for i = 1:length(Params.aggnames)
        eval(['xaggstst(' num2str(i) ') = ' Params.aggnames{i} ';']);
    end
    
    welfare = u * log(Params.b) - log(1-u+u*Params.b) ...
    + log(C) - (1-u)*h^(1+Params.gamma)/(1+Params.gamma) ...
    - (u+Params.upsilon*(1-u)) * q^(1+Params.kappa)/(1+Params.kappa) ...
    + Params.chi * log(G) - Params.psy*u;  

    welfare = welfare / (1-Params.beta);
end

