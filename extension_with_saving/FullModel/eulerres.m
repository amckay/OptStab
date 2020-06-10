function res = eulerres(par,parnext,beta,R,Rnext,wage,wagenext,dividend,dividendnext,lambda,lambdanext,M,Mnext,xthis)

global Params;


npp = Params.npp;

par = par2wide(par);
parnext = par2wide(parnext);


%figure out if we are doing a test of residuals off the grid
IsATest = exist('xthis','var');

if ~IsATest
    nc = Params.nc;
    nv = Params.nv;
    nn = Params.nn;
else
    nc = length(xthis);
    nv = nc;
    nn = nc;
end





for ip=1:npp
    
    N = nspline(par(Params.par_nind,Params.Nind(ip)));
    Q = nspline(par(Params.par_nind,Params.Qind(ip)));
    
    S = savingspline(par(Params.par_sind,ip));
    V = vspline(par(Params.par_vind,ip));
    
    
    if IsATest % xthis given
        sthis = interp_savspline(S,xthis);
    else % take x as starting value of assets
        xthis = S.x(1:end-1);  % last point in S is for interpolation
        sthis = S.y(1:end-1);
    end
    
    
    
    %we have end of period assets in last period and this period,
    %and we have to figure out this period's consumption:
    cthis = get_c(xthis,sthis,N,R,wage,dividend,lambda,ip);
    
    
    if(any(cthis<0))  %signal inadmissible value to routine 'broydn';
        %         disp('A');
        %         ip
        %         I = find(cthis<0)
        %         nthis(I)
        %         xthis(I)
        %         sthis(I)
        res = 1e100;
        return;
    end
    
    
    if(ip==1) % initialize z:
        res = initsize(zeros(nc+nn+nv,Params.npp),par,parnext,R,Rnext,wage,wagenext,dividend,dividendnext,lambda,lambdanext,M,Mnext,xthis);
    end
    
    
    assets = sthis;
    MUexp = 0;
    for jp = 1:npp 
        Sn  = savingspline(parnext(Params.par_sind,jp));
        Nn = nspline(parnext(Params.par_nind,Params.Nind(jp)));
        Qn = nspline(parnext(Params.par_nind,Params.Qind(jp)));

        qnext = interp_nspline(Qn,assets,true); %true -> take max of zero
        pp = transProb(ip,jp,Mnext,qnext,true);
    
        savingsnext = interp_savspline(Sn,assets);
        cnext = get_c(assets,savingsnext,Nn,Rnext,wagenext,dividendnext,lambdanext,jp);
    
        MUexp = MUexp + pp.*1./cnext;
    end
    
    
    
    
    % Euler residual:
    res(nv+nn+1:end,ip) = 1 - cthis.*(beta(ip)*Rnext*MUexp);
    
    clear nthis sthis cthis Sn Nn Qn qnext savingsnext cnext ;
    
    %--calculations for labor supply / search residuals--

    if Params.employed(ip)
        
        if IsATest % xthis given
            nthis = interp_nspline(N,xthis,false);
        else % take x as starting value of assets
            xthis = N.x(1:end-1);  % last point in N is for interpolation
            nthis = N.y(1:end-1);
        end
        sthis = interp_savspline(S,xthis);
        
        cthis = get_c(xthis,sthis,N,R,wage,dividend,lambda,ip);
        
        if(any(cthis<0))  %signal inadmissible value to routine 'broydn';
            %           disp('C');
            res = 1e100;
            return;
        end
        
        laborResid = (1-Params.tau)*lambda * (wage * nthis + dividend).^(-Params.tau) * wage*Params.Skill(ip).^(1-Params.tau) - cthis.* nthis.^Params.gamma;

    else
        
        if IsATest % xthis given
            qthis = interp_nspline(Q,xthis,false);
        else % take x as starting value of assets
            xthis = Q.x(1:end-1);  % last point in N is for interpolation
            qthis = Q.y(1:end-1);
        end
        
        V1 = vspline(par(1:Params.nv,Params.Nind(ip)));
        V2 = vspline(par(1:Params.nv,Params.Qind(ip)));
        dValue = interp_vspline(V1,xthis) - interp_vspline(V2,xthis);
        
        
        laborResid = (M*dValue).^(1./Params.kappa)./qthis -1;
        
        
    end
    
    %store results
    res(nv+1:nv+nn,ip) = laborResid;
    
    
    clear nthis1 sthis cthis;
    
    
    % ---------- calculations for V residuals --------------

    if IsATest % xthis given
        Vthis = interp_vspline(V,xthis);
    else
        Vthis = V.y(1:end-1);
        xthis_chv = V.x(1:end-1);
        xthis = exp(xthis_chv)-1;  %change of variables
    end
    sthis = interp_savspline(S,xthis);
    if any(sthis < Params.bhmin)
        res = 1e100;
        return;
    end
    
    if Params.employed(ip)
        nthis = interp_nspline(nspline(par(Params.par_nind,ip)),xthis,true);
    else
        nthis = 0;
    end

    cthis = get_c(xthis,sthis,N,R,wage,dividend,lambda,ip);
    assert(all(cthis > 0), 'error in eulerres: cthis is zero or negative.')


    %build the expected continuatiuon value
    assets = sthis;
    Vexp = 0;
    for jp = 1:npp
        if Params.employed(jp)
            Vcon =  interp_vspline(vspline(par(Params.par_vind,Params.Nind(jp))),assets);
        else
            V0 = interp_vspline(vspline(par(Params.par_vind,Params.Qind(jp))),assets);
            Vcon = V0 + Params.kappa/(1+Params.kappa) ...
                    * (interp_vspline(vspline(par(Params.par_vind,Params.Nind(jp))),assets) ...
                        -V0 ).^(1+1/Params.kappa) ;
            
        end
        
        pp = transProb(ip,jp,Mnext,0.0);
        Vexp = Vexp + pp .* Vcon;
    end
    
    
    res(1:nv,ip) = -Vthis +  log(cthis) - nthis.^(1+Params.gamma)/(1+Params.gamma) - Params.psy*~Params.employed(ip) + beta(ip)*Vexp;
    
    
end  %close loop over this period's income state



res = [reshape(res(1:nv,:),nv*Params.npp,1);
        reshape(res(nv+1:nv+nn,:),nn*Params.npp,1);
        reshape(res(nv+nn+1:end,:),nc*Params.npp,1)];

try
if max(abs(imag(res(:)))) > 1e-10
    res = 1e100;
    return;
end
end


end  % end function
