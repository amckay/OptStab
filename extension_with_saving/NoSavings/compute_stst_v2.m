function  [res, xaggstst, Env, X, lambda, welfare] = compute_stst_v2(X,Mode,varargin)
%Finds the steady state of the model
% MODE = 2 -> w is fixed and R and u are endogenous
% Remaining arguments are initial guesses / targets

global Params;


if strcmp(Mode,'calibration')

    lambda = X(1);
    h = X(2); % average hours among employed
    q = X(3); % average search effort among unemployed
    Vn = X(4);
    M = X(5);

    utarget = varargin{1};
    
    wage = 1/Params.mu  - Params.upsilon * Params.psi_1 *M^Params.psi_2 / h;

else

    lambda = X(1);
    h = X(2); % average hours among employed
    q = X(3); % average search effort among unemployed
    Vn = X(4);
    wage = Params.wbar;
    M = ((1/Params.mu  - wage)*h/Params.psi_1/Params.upsilon)^(1/Params.psi_2);
    

end

[u, ~, Y, dividend,J] = compute_steady_state_supply_side(h,q,wage);

    
p11 = 1-Params.upsilon + Params.upsilon * q * M;
p01 = Params.upsilon * (1- q *M);
% 1/ c_1 = beta * R * ( p11/c_1 + p01 / c_0)
% 1 = beta * R * ( p11 c_1/c_1 + p01 c_1 / c_0)
% 1 = beta * R * ( p11  + p01 / b )
R = 1 / Params.beta /(p11 + p01/Params.b);


[ cE, cU, h2, q2, Vn2 ] = SolveHousehold( h, Vn, lambda, dividend, M, wage );


% government budget constraint
C = (1-u)*cE + u*cU;
G = Params.chi * C;
lambda2 = (Y-J - G - (R-1)*Params.Bbar)/(C/lambda);


res = [lambda2-lambda; h2-h; q2- q; Vn2-Vn];

if strcmp(Mode,'calibration')
    
    res = [res;u-utarget];
end




if nargout >1

    % aggregate resource constraint (should clear by Walras)
    testResCon = Y-J-C-G;
    assert(abs(testResCon) < 1e-4*Y)


    [xaggstst, welfare] = agg_stst(R,lambda,h,q,Vn,wage);


    Env.R = R;
    Env.aggstst = xaggstst;



    X = xaggstst;


end

end % function

