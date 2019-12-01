function p=normmix_posterior(x,theta,prior,gt)
% NORMMIX_POSTERIOR Compute class probs. for a Gaussian mixture model.
%
% p=normmix_em(x,theta,prior)
% p=normmix_em(x,theta,prior,gt)
%
% x: n-by-d matrix
% theta{k}.mu: 1-by-d matrix, class expected value.
% theta{k}.Sigma: d-by-d matrix, class covariance.
% prior: 1-by-K matrix, the class probabilities
% gt: n-by-1 matrix with ground truth data (1..K indicating apriori known
%     class for a pixel, 0 indicating unknown class). NOT IMPLEMENTED!
%
% p: The posterior class probabilities, n-by-K matrix.

% $Id: normmix_posterior.m 4824 2014-12-07 21:56:53Z johanl $

if (nargin<4), gt = []; end

n = size(x,1);
d = size(x,2);
K = length(theta);

% Calculate p(\omega_i=k|x_i) \propto p(x_i|\omega_i=k) \pi_k
% for each pixel  i  and class  k
p = zeros(n,K); 
for k=1:K
  y = x-repmat(theta{k}.mu,[n,1]);
  p(:,k) = exp(-0.5*sum( ((y/theta{k}.Sigma).*y) ,2) ) / ...
	   ((2*pi)^(d/2)*det(theta{k}.Sigma)^0.5);
end
p = p*diag(prior);
p = p./repmat(sum(p,2),[1,K]);    

if (~isempty(gt))
  % Modify the posterior probabilities for pixels with known
  % "ground truth", gt.
  
  % Not implemented! (see Home Assignment)
  warning('fmsn20:gtNotImplementedInNormmixPosterior',...
          ['Ground truth available, but no implementation is ',...
           'available in normmix_posterior!'])
end
