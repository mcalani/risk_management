function nll = nll_lin( theta, data , vec)
% input:    (i) theta, coef vector, last element is error sd
%           (ii), data matrix, col1: y cols2:end: explan. variables (incl constant)
%           (iii), 0 = if vector of loglikelihoods, and 1 if sum should be
%           returned
beta = theta(1:end-1);
sig  = abs(theta(end))+0.000001;    % this ensures a non-zero variance
y = data(:,1);
x = data(:,2:end);
res = y - x*beta;
nll = (-0.5)*(-log(2*pi)-log(sig^2)-(res.^2/(sig^2)));

if vec
    nll = sum(nll);
end

end