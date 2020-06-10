function [u, M, Y, dividend,J] = compute_steady_state_supply_side(h,q,w)
% inputs average hours and average search effort

global Params;



A = 1;
% w = 1/Params.mu  - Params.upsilon * Params.psi_1 *M^Params.psi_2 / h;
M = ((1/Params.mu  - w)*h/Params.psi_1/Params.upsilon)^(1/Params.psi_2);
u = Params.upsilon*(1 - q * M)/(Params.upsilon*(1 - q * M) + q*M);
Hires = Params.upsilon.*(1-u);
Y = A * h * (1-u);
J = Params.psi_1 * M.^(Params.psi_2) * Hires;
dividend = (Y-J)/(1-u) - w * h;


end



