% calculations invariant distribution from transition matrix Pi
% Pi(i,j) is probability to go to i from j
function D = invdistr(Pi)
  assert(all(abs(sum(Pi)-1)<1e-10));
  opts.disp=0;
  [x,eval] = eigs(Pi,[],1,1+1e-10,opts);
  % [x,eval] = eigs(Pi,[],1,'lm',opts);
  assert(abs(eval-1)<1e-10);
  D = x/sum(x);
  assert(min(D)>-1e-12);
  D = max(D,0);
  D = D/sum(D);
