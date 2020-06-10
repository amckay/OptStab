%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION THAT DEFINES EULER RESIDUALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   par:           parameter vector for_ current savings function_
%   parnext:       parameter vector for_ next periods' savings function_
%   R,wage,Transf: 1+interest,wage,transfers
%   xthis:        grid at which_ to compute Euler residuals; default: knot points of spline
% output:
function z = eulerres(par,parnext,R,wage,Transf,xthis)
  global Params;

  npp = Params.npp;
  par = reshape(par,Params.nc,npp);
  parnext = reshape(parnext,Params.nc,npp);
  for ip=1:npp
    S = savingspline(par(:,ip));
    if(any(S.y(:)>S.x(:)))  %signal inadmissible value to routine 'broydn';
      z = 1e100;
      return;
    end

    if(nargin>5) % xthis given
      sthis = interp_savspline(S,xthis);
    else % take x as starting value for_ cash on hand
      xthis = S.x(1:end-1);  % last point in S is for interpolation
      sthis = S.y(1:end-1);
    end
    C = xthis-sthis;
    if(ip==1) % initialize z:
      z = initsize(zeros(length(C),Params.npp),par,parnext,R,wage,Transf,xthis);
    end

    assets = R*sthis + Transf;
    MU = margutil(C);
    MUexp = 0;

    for jp=1:npp
      pp = Params.transpp(ip,jp);
      if(pp>0)
	xnext = repmat(assets,1,size(Params.xY,1)) + wage*Params.ypp(jp)*repmat(Params.xY',size(xthis,1),1);
	Sn = savingspline(parnext(:,jp));
	if(any(Sn.y(:)>Sn.x(:)))
	  z = 1e100;
	  return;
	end
	ss = interp_savspline(Sn,xnext);
	Cnext = xnext - ss;
	
	MUnext = margutil(Cnext);
	MUexp = MUexp + pp*(MUnext*Params.wY);
      end
    end

    % Euler residual:
    % expressed in relative consumption units:
    z(:,ip) = 1 - invmargutil(Params.beta*R*MUexp)./C;
  end
  z = z(:);

