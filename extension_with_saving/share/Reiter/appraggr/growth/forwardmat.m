% compute transition matrix
% Pi(i,j) is probability to go from state j to state i!
function Pi = forwardmat(iFullMat,par,Transf,R,wage);
  global Params;

  nd = Params.ndstst+1;
  IR = []; IC = []; VV = [];
  for(jp=1:Params.npp)  % next period's state
    ir{jp} = [];
    ic{jp} = [];
    vv{jp} = [];
    S = savingspline(par(:,jp));
    for i=1:length(Params.xY)
      X = R*Params.knotDistrK + wage*Params.ypp(jp)*Params.xY(i) + Transf;
      Kend = interp_savspline(S,X);
      Ti = lineartrans(Params.knotDistrK,Kend);
      ir{jp} = [ir{jp};Ti.iTo];  % row indicates to which position we go
      ic{jp} = [ic{jp};Ti.iFr];  % column indicates which position we come from
      vv{jp} = [vv{jp};Params.wY(i)*Ti.Val];
    end
    if(iFullMat)
      offsj = (jp-1)*nd;
      for(ip=1:Params.npp)  % this period's state
	pp = Params.transpp(ip,jp);
	if(pp>0)
	  offsi = (ip-1)*nd;
	  IR = [IR;offsj+ir{jp}];  % where to go! take offsj
	  IC = [IC;offsi+ic{jp}];  % where come from! take offsi
	  VV = [VV;pp*vv{jp}];  % where come from! take offsi
	end
      end
    else
      Pij{jp} = sparse(ir{jp},ic{jp},vv{jp},nd,nd);
    end
  end
  if(iFullMat)
    nn = nd*Params.npp;
    Pi = sparse(IR,IC,VV,nn,nn);
  else
    Pi = op_concat(blockdiag(Pij),op_kron(Params.transpp',speye(nd)));
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION TO COMPUTE TRANSITION PROBABILITY MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = lineartrans(kgrid,k0)
  n = length(kgrid);
  k0 = max(min(k0,kgrid(n)),kgrid(1));
  iPos = lookup(kgrid,k0,0);
  iPos = reshape(iPos,1,length(k0));  % make sure it is column vector
  iPos = min(iPos,n-1);
  pHigh = (k0-kgrid(iPos))./(kgrid(iPos+1)-kgrid(iPos));
  pHigh = reshape(pHigh,1,length(k0));  % make sure it is column vector
  S.iFr = reshape(repmat(1:n,2,1),2*n,1);
  S.iTo = reshape([iPos;iPos+1],2*n,1);
  S.Val = reshape([1-pHigh;pHigh],2*n,1);

