function display(S)
  disp(sprintf('Size: (%d,%d)',S.n(1),S.n(2)));
  [ir,ic] = ind2sub2(S.n,S.i);
  disp([ir ic S.v]);
