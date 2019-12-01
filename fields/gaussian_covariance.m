function r=gaussian_covariance(dist,sigma2,rho,I,J)
% gaussian_covariance  calculates gaussian covariances
%
% r=gaussian_covariance(dist,sigma2,rho)
%
% r = matrix of covariances, calculated for  the distances in  dist
%     r has the same shape as dist
% sigma2, rho = the gaussian covariance parameters
%
% Approximate range = rho

% $Id: gaussian_covariance.m 5088 2017-11-04 11:03:21Z johanl $

if(nargin < 4), I=[]; J=[]; end

absrho = abs(rho);
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
%scale the distance
dunique = (dunique/absrho).^2;
%compute covariance
B = sigma2*exp(-2*dunique);
r(dpos) = B(J);
