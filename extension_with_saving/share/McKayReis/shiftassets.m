function Dshifted = scaleassets( D, shift )
% Takes asset distribution D over Params.knotDistrK and shifts asset
% holdings by shift.  The new distribution is 
% y = x + shift
% where x is distributed according to D.
% Just shift those with needy status

global Params;





  nd = Params.ndstst;
  IR = []; IC = []; VV = [];
  for jp=1:Params.npp  %just give to the needy
      
     
    if jp > Params.npp*2/3
        Kend = Params.knotDistrK + shift;
    else
        Kend = Params.knotDistrK;
    end
        
   
    
    Ti = lineartrans(Params.knotDistrK,Kend);
    ir{jp} = Ti.iTo;  % row indicates to which position we go
    ic{jp} = Ti.iFr;  % column indicates which position we come from
    vv{jp} = Ti.Val;
    
    
    Pij{jp} = sparse(ir{jp},ic{jp},vv{jp},nd,nd);
    
  end
  

Pi = blockdiag(Pij);

Dshifted = Pi * D;
  

end