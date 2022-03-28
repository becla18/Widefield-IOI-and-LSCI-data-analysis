function [res,gs] = gsr(Y,cst)
%GSR performs global signal regression on a time series matrix.
%   ARGUMENTS:
%   Y: data matrix, of size (n x d), where n is the no. of time points and d the number of
%   observations.
%   cst: optional, specify wether to include a constant term in the regression (default 1)
%   OUTPUTS:
%   res: residues of the regression, corresponding to the data matrix with the global signal
%   regressed out.
%   gs: Global signal vector
gs = mean(Y,2);
if nargin == 1
    cst = 1;
end
if cst
    X = [ones(size(gs)) gs];
else
    X = gs;
end
invX = pinv(X);
beta = invX*Y;
res = Y-X*beta;
end

