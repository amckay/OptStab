function Dnext = forward(D,Pi);
  % Dnext = D + Pi*D;
  Dnext = Pi*D;  % don't subtract eye in computing Pi!
