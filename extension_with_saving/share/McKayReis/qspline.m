function ns = qspline( par )

 global Params;

 
  
 
  ns.x = Params.ngrid;
  ns.y = par;
  
  
    
  
  %add end point 1e8
  ns.x = [ns.x; 1e8];
  ns.y = [ns.y; 0.5];
  

  ns.Slope = diff(ns.y) ./ diff(ns.x);
  
  

end

