function alpha_post=mrf_gaussian_post(alpha,theta,y)
% MRF_GAUSSIAN_POST Posterior alpha-parameters in a MRF with Gaussian data
%
%  alpha_post=mrf_gaussian_post(alpha,theta,y)
%
%  alpha : 1x1 or 1xK or mxnxK, the prior alpha-parameters.
%          See mrf_sim for the model specification.
%  theta : cell array of length K specifying the expectation vector and
%          covariance matrix for each k=1,...,K, as
%          theta{k}.mu and theta{k}.Sigma
%  y     : mxnxd, data image with d-dimensional Gaussian vector data.
%  alpha_post : mxnxK, the posterior alpha-parameters.

% Formula:
%   alpha_post_ik = alpha_ik + log(p(y_i|theta_k))

% $Id: mrf_gaussian_post.m 4837 2014-12-10 11:13:39Z johanl $

sz = size(y);
K = length(theta);
sz = [sz(1:2),K];
alpha = make_K_im(alpha,sz);
yc = colstack(y);
alpha_post = colstack(alpha);
for k=1:K
  ycm = yc-ones(sz(1)*sz(2),1)*theta{k}.mu(:)';
  alpha_post(:,k) = alpha_post(:,k) ...
      -1/2*log(2*pi*det(theta{k}.Sigma)) ...
      -0.5*sum(ycm.*(ycm/theta{k}.Sigma),2);
end
alpha_post = icolstack(alpha_post,sz(1:2));


function im=make_K_im(v,sz)
if (length(v)==1)
  im = v*ones(sz);
elseif ((size(v,1)==sz(1)) && (size(v,2)==sz(2)))
  if (size(v,3)==1)
    im = repmat(reshape(v,[sz(1:2),1]),[1,1,sz(3)]);
  else % v should already be the correct size; do nothing.
    im = v;
  end
else
  im = repmat(reshape(v,[1,1,sz(3)]),[sz(1:2),1]);
end
