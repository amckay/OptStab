
function xOut=plus(x1,x2)
  n1 = size(x1);
  n2 = size(x2);
  if(prod(n1)==1 & prod(n2)>1)
      x1 = repmat(x1,n2(1),n2(2));
  end
   if(prod(n2)==1 & prod(n1)>1)
      x2 = repmat(x2,n1(1),n1(2));
  end
  x = getval(x1);
  y = getval(x2);
  dfdx1 = 1;
  dfdx2 = 1;
  dfdx1 = dfdx1(:);
  dfdx2 = dfdx2(:);
  if(notisa(x1,'deriv1'))
    np = nindep(x2);
  else
    np = nindep(x1);
  end
  if(notisa(x1,'deriv1'))
    xOut.v = plus(x1,x2.v);
    xOut.d = elelmult_eachcol(dfdx2,x2.d);
  elseif(notisa(x2,'deriv1'))
    xOut.v = plus(x1.v,x2);
    xOut.d = elelmult_eachcol(dfdx1,x1.d);
  else  % both deriv1
    n1d = size(x1.d);
    np = n1d(end);
    xOut.v = plus(x1.v,x2.v);
    xOut.d = elelmult_eachcol(dfdx1,x1.d) + elelmult_eachcol(dfdx2,x2.d);
  end
  xOut=class(xOut,'deriv1');
