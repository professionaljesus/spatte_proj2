function r=matern_covariance(dist,sigma2,kappa,nu,I,J)
% matern_covariance  calculates Matérn covariances
%
% r=matern_covariance(dist,sigma2,kappa,nu)
%
% r = matrix of covariances, calculated for  the distances in  dist
%     r has the same shape as dist
% sigma2, kappa, nu = the Matérn covariance parameters
%
% Approximate range = sqrt(8*nu)/kappa

% $Id: matern_covariance.m 5088 2017-11-04 11:03:21Z johanl $

if(nargin < 5), I=[]; J=[]; end

absnu = abs(nu);
absk = abs(kappa);
%compute range
rho = sqrt(8*absnu)/kappa;
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
B = log(besselk(absnu,absk*dunique));
B = exp(log(sigma2)-gammaln(absnu)-(absnu-1)*log(2) + ...
    absnu.*log(absk*dunique) + B);
%sanity check, large nu values will need replacement for small dunique
%values
if any(isinf(B))
    B(isinf(B)) = sigma2*exp(-2*(dunique(isinf(B))/rho).^2);
end
r(dpos) = B(J);
