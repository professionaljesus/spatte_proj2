function [mu, Sigma] = gibbs_mu_sigma(y)
% GIBBS_MU_SIGMA Sample the posterior mean and covariance given y
%
% [mu, Sigma] = gibbs_mu_sigma(y)
%
% x: n-by-d matrix
%
% mu: 1-by-d sample from the posterior mean
% Sigma: d-by-d sample from the posterior covariance

% $Id: gibbs_mu_sigma.m 4824 2014-12-07 21:56:53Z johanl $

%compute mean and Syy
[n,d] = size(y);
mu = mean(y);
y_m =  bsxfun(@minus, y, mu);
Syy = y_m'*y_m;

%sample posterior Sigma (ensure valid wishart for small values of n)
%TODO: USE A PROPER PRIOR!
if n<=d, Syy = Syy + 1e-5*eye(size(Syy)); end
Sigma = iwishrnd(Syy, max(n-d-2,d+1) );
%sample posterior mean
R = chol(Sigma/n);
mu = mu + reshape(R'*randn(d,1), size(mu)); 
