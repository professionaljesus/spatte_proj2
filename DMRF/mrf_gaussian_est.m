function theta_hat=mrf_gaussian_est(W,y)
% MRF_GAUSSIAN_EST Weighted estimate of Gaussian parameters
%
%  theta_est=mrf_gaussian_est(W,y)
%
%  W : mxnxK, weights for each data point.
%      one set of parameters is estimated for every k=1,...,K
%  y : mxnxd, data image with d-dimensional Gaussian vector data.
%  theta_hat : cell array of length K with the estimated
%              expectation vector and covariance matrix for
%              each k=1,...,K, as theta{k}.mu and theta{k}.Sigma

% $Id: mrf_gaussian_est.m 4837 2014-12-10 11:13:39Z johanl $

szW = size(W);
szy = size(y);
if ((szW(1)~=szy(1)) || (szW(2)~=szy(2)))
  error('The first two dimensions of W and y must have the same size.')
end
K = size(W,3);
d = size(y,3);
theta_hat = cell(1,K);
Wc = colstack(W);
yc = colstack(y);
for k=1:K
  sW = sum(Wc(:,k));
  theta_hat{k}.mu = (Wc(:,k)'*yc)/sW;
  yc_ = yc-ones(szy(1)*szy(2),1)*theta_hat{k}.mu;
  theta_hat{k}.Sigma = (yc_'*((Wc(:,k)*ones(1,d)).*yc_))/sW;
end
