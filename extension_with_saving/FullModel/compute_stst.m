function  [xaggstst, par, D, Env, X] = compute_stst(MODE,ha,qa,lambdaa,Ra,Ba,wa)
%Finds the steady state of the model
% MODE = 1 -> R and w are fixed and B and u are endogenous
% MODE = 2 -> B and w are fixed and R and u are endogenous
% Remaining arguments are initial guesses / targets


global Params;

MYTOL = 1e-6;


if MODE == 1
    R = 1.005;
    wage = wa;
elseif MODE == 2
    Btarg = Ba;
    wage = wa;
    R = Ra; % initial guess
    
else
    error('Mode not recognized')
end


lambda = lambdaa;
h = ha; % average hours among employed
q = qa; % average search effort among unemployed



if MODE ==1
    for itlam = 1:100
        for ithq = 1:100
            % given h, q and u we can figure out M and then from the firm's FOC we
            % can figure out w
%             M = (1-u/Params.upsilon)/q ;
%             wage =   1/Params.mu  - Params.psi_1*M^Params.psi_2 / h;
% 
%             [uimp, M, Y, dividend,J] = compute_steady_state_supply_side(h,q,wage);
%             assert(abs(u - uimp) < 1e-5)
            
            [u, M, Y, dividend,J] = compute_steady_state_supply_side(h,q,wage);

            % SOLVE FOR THE CONSUMPTION FUNCTION BY COLLOCATION:
            % high precision in solution of savings function_ is required to achieve
            %   converge in outer loop!
            [par,check] = broydn(@eulerres_stst,Params.parstart,[1e-11,0,1],R,wage,dividend,lambda,M);
            if(check~=0)
                %save broyError par;
                warning('compute_stst:broyerror','broydn not converged');
            end


            Params.parstart = par;
            par = par2wide(par);

            Pi = forwardmat(par,M);
            D = invdistr(Pi);
            % D = invd2(Pi);
            clear Pi;

            % loop to convergence on h and q, checking consistency
            hq = expect_L(D,par,M);

            testhq = max(abs(hq - [h q]));
            if testhq < MYTOL, break; end
            h = hq(1);
            q = hq(2);

        end %hq
        
        % government budget constraint
        ET = expect_tax(D,par,wage,dividend,M);
        [C, EInc] = expect_C(D,par,R,wage,dividend,lambda,M);
        B = expect_k(D);
        assert(abs(C - EInc - (R-1)*B) < 50*MYTOL)
        G = Params.chi * C;
%         if strcmp(ADJ,'R')
%           lambda2 = (ET(3) - G - (R-1)*Params.Bbar)/(ET(1) + Params.b*ET(2));
%         else
          lambda2 = (ET(3) - G - (R-1)*B)/(ET(1) + Params.b*ET(2));
%         end
        

        testlam = abs(lambda2 - lambda);
        if testlam < MYTOL, break; end
        lambda = lambda2;

    end % lambda
    assert(abs(EInc - lambda * (ET(1) + Params.b*ET(2))) < 1e-8)

elseif MODE == 2
            

    for itR = 1:100
        for itlam = 1:100
            for ithq = 1:100
                [u, M, Y, dividend,J] = compute_steady_state_supply_side(h,q);


                % GIVEN K, SOLVE FOR THE CONSUMPTION FUNCTION BY COLLOCATION:
                % high precision in solution of savings function_ is required to achieve
                %   converge in outer loop!
                [par,check] = broydn(@eulerres_stst,Params.parstart,[1e-11,0,1],R,wage,dividend,lambda,M);
                if(check~=0)
                    %save broyError par;
                    warning('compute_stst:broyerror','broydn not converged');
                end


                Params.parstart = par;
                par = par2wide(par);

                Pi = forwardmat(par,M);
                D = invdistr(Pi);
                % D = invd2(Pi);
                clear Pi;

                % loop to convergence on h and q, checking consistency
                hq = expect_L(D,par,M);

                testhq = max(abs(hq - [h q]));
                if testhq < MYTOL, break; end
                h = hq(1);
                q = hq(2);

            end %hq

            % government budget constraint
            ET = expect_tax(D,par,wage,dividend,M);
            [C, EInc] = expect_C(D,par,R,wage,dividend,lambda,M);
            B = expect_k(D);
            assert(abs(C - EInc - (R-1)*B) < 50*MYTOL)
            G = Params.chi * C;
    %         if strcmp(ADJ,'R')
    %           lambda2 = (ET(3) - G - (R-1)*Params.Bbar)/(ET(1) + Params.b*ET(2));
    %         else
              lambda2 = (ET(3) - G - (R-1)*B)/(ET(1) + Params.b*ET(2));
    %         end


            testlam = abs(lambda2 - lambda);
            if testlam < MYTOL, break; end
            lambda = lambda2;

        end % lambda
        assert(abs(EInc - lambda * (ET(1) + Params.b*ET(2))) < 1e-8)

        % bond market clearing
        B = expect_k(D);
        testB = Btarg-B
        if abs(testB) < MYTOL, break; end
        
        if abs(testB) > 0.05    
            R = R + 0.0005*testB;
        else
            R = R + 0.05*testB;
        end
        bps = 1e4*(R-1)
        
    end % it over R

end

% aggregate resource constraint (should clear by Walras)
testResCon = Y-J-C-G;
assert(abs(testResCon) < 100*MYTOL)


xaggstst= agg_stst(wage,R,lambda,h,q,D,par);

if(nargout>=4)
    Env.Dstst = D;
    Env.R = R;
    Env.wage = wage;
    Env.parstst = par;
    Env.aggstst = xaggstst;
end

if (nargout >= 5)
    X = [zeros(Params.nStatesDistr-1,1);xaggstst;par2long(par)];
end



end % function
