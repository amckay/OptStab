function vs = vspline( par )

 global Params;

 if isempty(par)
     vs = [];
     return;
 end
  
 
  vs.x = Params.vgrid;
  vs.y = par;
  
  
    
  
  %add end point 1e8
  vs.x = [vs.x; 1e8];
  vs.y = [vs.y; 0];
  

  vs.Slope = diff(vs.y) ./ diff(vs.x);
  
  

end

