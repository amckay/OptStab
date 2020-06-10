function xaggstst= agg_stst(C, u, R, Y, G, mathcalH, wage, M, J, B, dividend, lambda) %#ok


    global Params;
   
    
    A = 1; %#ok
    Gshock = 1; %#ok
    mpshock = 1;%#ok
    
    I = R; 
    ppi = 1; %#ok
    pstar = 1; %#ok

    
    pbarB = Y/(1-  (1-Params.theta)/I);
    pbarA = pbarB; %#ok
    
    xaggstst = zeros(Params.nagg,1);
    for i = 1:length(Params.aggnames)
        eval(['xaggstst(' num2str(i) ') = ' Params.aggnames{i} ';']);
    end
    
    
end

