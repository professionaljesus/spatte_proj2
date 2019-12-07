function [logp, D_logp, D2_logp]= proj2_gmrf_taylor_skeleton(x_0, y, A, Q, E, pars)
% GMRF_TAYLOR_SKELETON  Taylor expansion of the conditional for non-Gaussian observations
%
% [logp, D_logp, D2_logp]= GMRF_taylor_Po(x_0, y, A, Q, E)
%
% x_0 = value at which to compute taylor expansion, as a column with N elements
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-N
% Q = the precision matrix, sparse N-by-N
% E = the population count in each region, as a column with n elements
% pars = possibly additional parameters needed by the observation likelihood
%
% Function should return taylor expansion, gradient and Hessian.
%
% This is only a skeleton for Home Assignment 2.

% $Id: gmrf_taylor_skeleton.m 5107 2017-11-12 13:35:17Z johanl $

%ensure that E=1 if E not given (i.e. same/no population weight in all regions)
if nargin<5, E=1; end
if nargin<6, par=[]; end

%compute observed part of the field
z = A*x_0;
senior_x = = 
%compute log( p(y|z,theta) )
f = -log(factorial(y)) + y.*(log(E + z)) - E.*exp(z); %Poisson?

nabla_f = y - E.*exp(z);
Hessian_f = spdiag(E.*exp(z));
%compute -log p(x|y,theta)
logp = -(sum(f) + (y - z)'*nabla_f + 0.5*(y - z)'*Hessian_f*(y - z));

if nargout>1
  %compute derivatives (if needed, i.e. nargout>1)
  df = nabla_f;
  D_logp = [] - A'*df;
end

if nargout>2
  %compute hessian (if needed, i.e. nargout>2)
  d2f = Hessian_f;
  n = size(A,1);
  D2_logp = -Q - A'*spdiags(d2f,0,n,n)*A;
end
