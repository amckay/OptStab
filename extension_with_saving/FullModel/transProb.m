function [prob ] = transProb(i,j,M,q,expect)

global Params;

if ~exist('expect','var')
    expect = false;
end



jf = M*q;


if Params.nbeta == 1
    js = Params.upsilon*(1-jf);
    if Params.employed(i)  % was employed
        if Params.employed(j) % still employed
            prob = 1-js;
        else  % lose job
            prob = js;
        end
    else  % was unemployed
        if Params.employed(j)  % found job
            prob = jf;
        else  % still unempl
            prob = 1-jf;
        end
    end
else
    js = Params.upsilon_groups(j)*(1-jf);
    if expect
        T1 = kron(Params.betaSwitchExpect, ones(2));
    else
        T1 = kron(Params.betaSwitch, ones(2));
    end
   
   if Params.employed(i)  % was employed
        if Params.employed(j) % still employed
            T2 = 1-js;
        else  % lose job
            T2 = js;
        end
    else  % was unemployed
        if Params.employed(j)  % found job
            T2 = jf;
        else  % still unempl
            T2 = 1-jf;
        end
    end

     prob = T1(i,j) * T2;
end
        
    
    
end

