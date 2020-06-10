% compute transition matrix
% Pi(i,j) is probability to go from state j to state i!
function Pi = forwardmat(par,M)



global Params;

nd = Params.ndstst;
offsi = ((1:Params.npp)-1)*nd;
offsj = ((1:Params.npp)-1)*nd;

IR = []; IC = []; VV = [];



for ip=1:Params.npp  % labor market state at end of previous period
    
    Q = nspline(par(Params.par_nind,Params.Qind(ip)));
    q = interp_nspline(Q,Params.knotDistrK,true);
    
    for jp = 1:Params.npp  % state after transitions
        pp = transProb(ip,jp,M,q);
   
        S = savingspline(par(Params.par_sind,jp));
        Kend = interp_savspline(S,Params.knotDistrK);
        
        
        Ti = lineartrans(Params.knotDistrK,Kend);
        IR = [IR; offsj(jp)+Ti.iTo];  % row indicates to which position we go to
        IC = [IC; offsi(ip)+Ti.iFr];  % column indicates which position we come from
        VV = [VV; reshape([pp';pp'],size(Ti.Val,1),1).*Ti.Val];
       
    end
        
    
end  %end loop over this period state


Pi = sparse(IR,IC,VV,nd * Params.npp,nd* Params.npp);


end