%{
AGM - Notes
X is the point at which to take derivatives. 
iDeriv is the index of which variables in X are to be differentiated with respect to.
Other parts of X are assumed constant.

nmax tells adjacob how many derivatives to take at a time.  I believe this has to do with
memory management.

I am not sure what the custom for whether X and func should be/produce row or column vectors.

%}

clc; clear all;
X = ones(1,2);
iDeriv = 1:2;
nmax = 1;

J = adjacob(@test_equ,X,iDeriv,nmax)