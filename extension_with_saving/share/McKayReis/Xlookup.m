function x = Xlookup(X,names, nm)
%look up index of X corresponding to names == nm

x = X(strcmp(names,nm));