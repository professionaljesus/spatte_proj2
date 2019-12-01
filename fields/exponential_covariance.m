function r=exponential_covariance(dist,sigma2,kappa,I,J)
% exponential_covariance  calculates exponential covariances
%
% r=exponential_covariance(dist,sigma2,kappa)
%
% r = matrix of covariances, calculated for  the distances in  dist
%     r has the same shape as dist
% sigma2, kappa = the exponential covariance parameters
%
% Approximate range = 2/kappa

% $Id: exponential_covariance.m 5088 2017-11-04 11:03:21Z johanl $

if(nargin < 4), I=[]; J=[]; end

absk = abs(kappa);
%compute range
dpos = (dist>0);
r = zeros(size(dist));
r(~dpos) = sigma2;
%% Less sensitive to large distances:
if( isempty(I) || isempty(J) )
  [dunique,~,J] = unique(dist(dpos));
else
  dunique = dist(dpos);
  dunique = dunique(I);
end
B = sigma2*exp(-absk*dunique);
r(dpos) = B(J);
