function [y,P,Pvar]=pca(x,use_eig,x_mean)
% PCA Principal component transformation of data.
%
%  y = pca(x)
%  [y,P,Pvar] = pca(x)
%
%  x: n d-dimenisonal data points as an n-by-d matrix.
%  y: Transformed data points.
%  P: Principal component axes, orthogonal d-by-d matrix.
%     Each column P(:,k) is the covariance eigenvector
%     corresponing to eigenvalue Pvar(k).
%  Pvar: Principal component variances, 1-by-d vector,
%        sorted by decreasing value.
%
%  y_i = (x_i-m)*P
%  Cov(x_i) = P*diag(Pvar)*P'
%  Cov(y_i) = diag(Pvar)
%
% See also: pca, markpca

% Copyright (c) 2001,2002 by Finn Lindgren
% $Revision: 3082 $  $Date: 2006-10-26 16:09:10 +0200 (tor, 26 okt 2006) $

% Predef:
%   m = mean(x,1)
%   z_i = x_i-m
%
% Theory, using diagonalisation of Sigma:
%   Sigma = z'*z/n;
%   [P,Pvar] = eig(Sigma); Pvar = diag(Pvar);
%   Cov(z_i) = Sigma = P*diag(Pvar)*P'
%   y = z*P  ==>  Cov(y_i) = diag(Pvar)
%   (Note that eig may not give the eigenvalues in the desired order)
%
% Another, using SVD factorisation directly on data:
%   [U,S,P] = svd(z,0);
%   z'*z = (U*S*P')'*U*S*P' = P*S'*U'*U*S*P' = P*S'*S*P'
%   y = z*P  ==>  y'*y = P'*z'*z*P = S'*S = n*diag(Pvar)
%   Cov(y_i) = diag(Pvar)

if (nargin<2)
  use_eig = [];
end
if isempty(use_eig)
  use_eig = 0;
end

if nargin<3
  m = mean(x,1);
  z = x-ones(size(x,1),1)*m;
else
  z = x-ones(size(x,1),1)*x_mean;
end

if (use_eig)
  Sigma = z'*z/size(x,1);
  [P,Pvar] = eig(Sigma);
  % Sort eigenvalues in decreasing order:
  [Pvar,idx] = sort(diag(Pvar));
  Pvar = Pvar(end:-1:1)';
  P = P(:,idx(end:-1:1));
else
  [U,S,P]=svd(z/sqrt(size(x,1)),0);
  Pvar = diag(S'*S);
end
y = z*P;
