function Env = dojacob(funcname,X,Env,njac,varargin)
  if(nargin<4)
    njac = 1000;
  end
  Env.jac0 = sparse(adjacob(funcname,X,Env.varix.x,njac,varargin{:}));
  Env.jac1 = sparse(adjacob(funcname,X,Env.varix.xlag,njac,varargin{:}));
  Env.jace = sparse(adjacob(funcname,X,Env.varix.eps,njac,varargin{:}));
  Env.jact = sparse(adjacob(funcname,X,Env.varix.eta,njac,varargin{:}));

%  assert(nnz(Env.jac1(Env.iEquBWS,Env.iVarDec))==0);
  if(isfield(Env,'iVarJump'))
    assert(nnz(Env.jac1(:,Env.iVarJump))==0);
  end
  



