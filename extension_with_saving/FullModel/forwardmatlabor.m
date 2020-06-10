% compute transition matrix
% Pi(i,j) is probability to go from state j to state i!
function Pi = forwardmatlabor(par,M)

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
   
        IR = [IR;offsj(jp)+(1:nd)'];  % where to go
        IC = [IC;offsi(ip)+(1:nd)'];  % where come from
        VV = [VV;pp];
       
    end
    
    
end  %end loop over this period state

Pi = sparse(IR,IC,VV,nd*Params.npp,nd*Params.npp);


end

