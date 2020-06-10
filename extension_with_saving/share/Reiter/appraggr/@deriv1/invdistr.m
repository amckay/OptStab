function D = invdistr(Pi)
  D.v = invdistr(Pi.v);

  n = length(Pi.v);
  p = nindep(Pi);
  for i=1:p
    rhs = reshape(Pi.d(:,i),n,n)*D.v;
  end

  D=class(D,'deriv1');
