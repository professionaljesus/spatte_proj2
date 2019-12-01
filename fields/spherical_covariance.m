function r=spherical_covariance(dist,sigma2,rho,I,J)
% spherical_covariance calculates spherical covariances
%
% r=spherical_covariance(dist,sigma2,rho)
%
% r = matrix of covariances, calculated for  the distances in  dist
%     r has the same shape as dist
% sigma2, rho = the spherical covariance parameters
%
% Approximate range = 0.73*rho

% $Id: spherical_covariance.m 5090 2017-11-05 19:35:05Z johanl $

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
dunique = dunique/absrho;
%compute covariance
B = sigma2 * (max(1-1.5*dunique+0.5*dunique.^3,0) .* (dunique<1));
r(dpos) = B(J);
