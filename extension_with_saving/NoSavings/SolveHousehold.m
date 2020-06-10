function [ cE, cU, h, q, Vn ] = SolveHousehold( h0, Vn0, lambda, dividend, M, wage )

    global Params;
    h = fzero(@SolveHouseholdInner_h,h0,[], lambda, dividend,wage);
    cE = lambda*(wage * h + dividend)^(1-Params.tau);
    cU = Params.b * cE;
    
    Vn = fzero(@SolveHouseholdInner_Vn,Vn0,[],h,M);
    
   
    q = (M .* Vn).^(1./Params.kappa);

end


function resid = SolveHouseholdInner_h(h, lambda, dividend,wage)

    global Params;

    cE = lambda*(wage * h + dividend)^(1-Params.tau);
    resid = lambda * (1-Params.tau) * (wage * h + dividend)^(-Params.tau) * wage - cE * h^Params.gamma;



end

function resid = SolveHouseholdInner_Vn(Vn, h, M)

    global Params;

    q = (M .* Vn).^(1./Params.kappa);
    
    Vn_impl = (-log(Params.b) - h.^(1+Params.gamma)/(1+Params.gamma)+Params.psy)/(1-Params.beta*(1-Params.upsilon)*(1-Params.kappa/(1+Params.kappa)*q*M));
    
    resid = Vn_impl - Vn;



end