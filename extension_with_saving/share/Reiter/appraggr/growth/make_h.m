function H = make_h(ii,n)
  global Params;
  if(nargin<2)
    n = Params.ndstst;
  end
  H = zeros(n,length(ii));
  for j=ii
    h = expect_k([],j);
    h = reshape(h,Params.ndstst+1,Params.npp);
    m = length(h)-1;
    h2 = h(2:m+1,:) - repmat(h(1,:),m,1);
    H(1:length(h2(:)),j) = h2(:);
  end
