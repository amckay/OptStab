function [parNew, Cprev, vnew] = egm(par,beta,R,wage,dividend,lambda,M,do_q)
%par = egm_c(par,R,Rnext,wage,wagenext,ppinext,labormarketstatus)
%iterate once on consumption and labor supply using the endog grid method

global Params;


npp = Params.npp;


%% Section 1: solve for savings rule using EGM

% asset grid
xthis = [0; Params.knotXi];


nassets = length(xthis);
nc = Params.nc;
assert(nc == nassets);

%compute consumption and marg util for each asset level and income level

MU = NaN(nassets,npp);
for jp=1:npp
    
    N = nspline(par(Params.par_nind,Params.Nind(jp)));
    S = savingspline(par(Params.par_sind,jp));
    sthis = interp_savspline(S,xthis);
    cthis = get_c(xthis,sthis,N,R,wage,dividend,lambda,jp);
    if ~ all(cthis>0)
        error('negative c')
    end
    assert(all(cthis>0));
    MU(:,jp) = 1./cthis;
end

%compute expected marg util
MUexp = zeros(nassets,npp);
for ip = 1:npp
    for jp = 1:npp    
        Q = qspline(par(Params.par_nind,Params.Qind(jp)));
        qthis = interp_nspline(Q,xthis,true); %true -> take max of zero
        pp = transProb(ip,jp,M,qthis,true);
        MUexp(:,ip) = MUexp(:,ip) + pp .* MU(:,jp);
    end
end
clear MU;

Cprev = 1./(bsxfun(@times,beta,R*MUexp));
assert(all(Cprev(:) > 0))


%at this stage we know the previous c and n on the b' grid, but we
%don't know n or b
assert(all(size(Cprev) == [Params.nc, npp]))
%bprev = egm_bn(repmat(xthis,1,Params.npp),nthis,Cprev,repmat(xthis,1,Params.npp),R-1,ppinext,wage);
bprev1 = repmat(xthis,1,npp);
Inc = zeros(nassets,npp);
for it = 1:100
    for ip = 1:npp
        Inc(:,ip) = get_laborIncome(bprev1(:,ip),N,wage,dividend,lambda,ip);
    end
    bprev = (repmat(xthis,1,Params.npp) + Cprev - Inc)/R; 
    if any(abs(imag(bprev(:))>1e-9))
        error('why imag')
    end
    test = max(abs(bprev(:) - bprev1(:)));
    if test < 1e-5, break; end
    bprev1 = bprev;
end


%pack results back into par
parNew =  NaN(size(par));

parNew(Params.par_sind(1),:) = bprev(1,:);
for ip = 1:npp
    
    XX = [bprev(:,ip); 1e8];
    tmp = (xthis(end)-xthis(end-1))/(bprev(end,ip)-bprev(end-1,ip))*(1e8-bprev(end,ip)) + xthis(end);
    YY = [xthis; tmp];
    
    parNew(Params.par_sind(2:end),ip) = interp1(XX,YY,bprev(1,ip) + Params.knotXi,'linear');
    
    %if ip == 3
    %    [bprev(end-8:end,ip)'; xthis(end-8:end)';bprev(1,ip) + Params.knotXi(end-8:end)'; parNew(Params.par_sind(end-8:end),ip)']
    %end
end

% 
% 
% % consistency check
% xthis_chv = Params.vgrid;
% xthis = exp(xthis_chv)-1;  %change of variables
% for ip = 1:npp
%     
%     S = savingspline(parNew(Params.par_sind,ip));
%     N = nspline(par(Params.par_nind,ip));
%     
%     sthis = interp_savspline(S,xthis);
%     [cthis, nthis, margTax, taxpaid, taxableIncome] = get_cnt(xthis,sthis*ppinext,N,R,wage,ip);
%     if ~all(cthis > 0)
%         [bprev(end,ip) xthis(end) bprev(1,ip)+Params.knotXi(end) parNew(Params.par_sind(end),ip)]
%         ip
%         I = cthis <= 0;
%         find(I)
%         xthis(I)
%         nthis(I)
%         sthis(I)
%         R
%         wage
%         taxableIncome(I)
%     end
%     assert(all(cthis > 0), 'error in egm: cthis is zero or negative.')
%     end

%% Section 2 --  Solve for n rule using static first order condition
% we have already done this before, but we may not have solved for the
% decision rule in all parts of the state space

xthis = Params.ngrid;

for ip = 1:npp
    if Params.employed(ip)
        S = savingspline(parNew(Params.par_sind,ip));
        savings = interp_savspline(S,xthis);

        N = nspline(par(Params.par_nind,ip));  %for an initial gues
        nnew = interp_nspline(N,xthis,true);
        nnew = solveN(nnew, xthis,savings,R,wage,dividend,lambda,ip);

        parNew(Params.par_nind,ip) = nnew;
    end
end

%% Section 3 -- Solve for V using Bellman equation
vnew = NaN(Params.nv,npp);


xthis_chv = Params.vgrid;
xthis = exp(xthis_chv)-1;  %change of variables


%loop over this period's income state
for ip = 1:npp

    N1 = nspline(parNew(Params.par_nind,Params.Nind(ip)));
    N1 = interp_nspline(N1,xthis,true);
    
    S = savingspline(parNew(Params.par_sind,ip));
    
    sthis = interp_savspline(S,xthis);
    cthis = get_c(xthis,sthis,N1,R,wage,dividend,lambda,ip);
    if ~all(cthis > 0)
        ip
        I = cthis <= 0;
        find(I)
        xthis(I)
        N(I)
        sthis(I)
        R
        wage
        lambda
        dividend
    end
    assert(all(cthis > 0), 'error in egm: cthis is zero or negative.')
    
    
    %build the expected continuatiuon value
    assets = sthis;
    Vexp = 0;
    Vcon = zeros(length(assets),npp);
    for jp = 1:npp
 
        if Params.employed(jp)
            Vcon(:,jp) =  interp_vspline(vspline(par(Params.par_vind,Params.Nind(jp))),assets);
        else
            V0 = interp_vspline(vspline(par(Params.par_vind,Params.Qind(jp))),assets);
            Vcon(:,jp) = V0 + Params.kappa/(1+Params.kappa) ...
                    * max(interp_vspline(vspline(par(Params.par_vind,Params.Nind(jp))),assets) ...
                        -V0,1e-2 ).^(1+1/Params.kappa) ;
            
        end
    end
    for jp = 1:npp
        
        pp = transProb(ip,jp,M,0.0);
        Vexp  = Vexp + pp .* Vcon(:,jp);
    end
    
    
    vnew(:,ip) = log(cthis) - Params.employed(ip)*(N1.^(1+Params.gamma))/(1+Params.gamma)...
        - ~Params.employed(ip)* Params.psy ...
        + beta(ip)*Vexp; % note I have left off public goods benefit since these are common to employed and unemployed
    
    
%     if ip ==2
%         save tmp.mat
%     end
end  %close loop over this period's income state


%pack results back into par
parNew(Params.par_vind,:) = vnew;


%% Section 4 Solve for q for searchers

if do_q
    for ip = 1:npp
        if ~Params.employed(ip)
            VE = interp_vspline(vspline(par(Params.par_vind,Params.Nind(ip))),Params.ngrid);
            VU = interp_vspline(vspline(par(Params.par_vind,Params.Qind(ip))),Params.ngrid);
            qnew = (M*(VE - VU)).^(1/Params.kappa);
            parNew(Params.par_nind,ip) = qnew;
        end
    end
end

end


function [n] = solveN(n0, b,savings,R,wage,dividend,lambda,incomeIndex)
% solves non-linear equation for labor given assets and savings


global Params;

n = n0;

for it = 1:100
    
    c = R * b  - savings + get_laborIncome(b,n,wage,dividend,lambda,incomeIndex);
    f2 = (1-Params.tau)*lambda * (wage * n + dividend).^(-Params.tau) * wage *Params.Skill(incomeIndex).^(1-Params.tau)...
     - c.* n.^Params.gamma;
    
    
    
    
    if all(abs(f2(:)) < 1e-14)
        break
    end
    
    J2n = -Params.tau * (1-Params.tau)*lambda * (wage * n + dividend).^(-Params.tau-1) * wage^2 *Params.Skill(incomeIndex).^(1-Params.tau) ...
        - Params.gamma * c.* n.^(Params.gamma-1) ...
        -n.^Params.gamma .* (1-Params.tau) .* lambda .* (wage * n + dividend).^(-Params.tau) * wage *Params.Skill(incomeIndex).^(1-Params.tau);
    
    
    D = -f2./J2n;
    
    n = n + D;
    
    
    
    
    
end

end

