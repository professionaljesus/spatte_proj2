function r=cauchy_covariance(dist,sigma2,rho,x,I,J)
% cauchy_covariance  calculates Matérn covariances
%
% r=cauchy_covariance(dist,sigma2,rho,x)
%
% r = matrix of covariances, calculated for  the distances in  dist
%     r has the same shape as dist
% sigma2, rho, x = the Matérn covariance parameters
%
% Approximate range = rho*sqrt(10^(1/x)-1)

% $Id: cauchy_covariance.m 5088 2017-11-04 11:03:21Z johanl $

if(nargin < 5), I=[]; J=[]; end

absrho = abs(rho);
n_absx = -abs(x);
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
%scale the distance
dunique = (dunique/absrho).^2;
%compute covariance
B = sigma2 * (1+dunique).^n_absx;
r(dpos) = B(J);

%x = d*sqrt( (1/10)^(1/x)-1 )