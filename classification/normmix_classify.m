function [cl,cl_ind,p]=normmix_classify(x,theta,prior,gt)
% NORMMIX_CLASSIFY Classify data in a Gaussian mixture model.
%
% [cl,cl_ind,p]=normmix_classify(x,theta,prior)
% [cl,cl_ind,p]=normmix_classify(x,theta,prior,gt)
% [cl,cl_ind]=normmix_classify(p)
%
% x: n-by-d matrix, data, one d-dimensional observation per row.
% theta{k}.mu: 1-by-d matrix
% theta{k}.Sigma: d-by-d matrix
% prior: 1-by-K matrix, the prior class probabilities
%
% cl: n-by-1 matrix, the classification indices
% cl_ind: n-by-K matrix, classification indicator matrix,
%         useful for classifications images with rgbimage
% p: n-by-K matrix, the posterior class probabilities, as computed
%    by normmix_posterior or normmix_em
% gt: ground truth data, pre-classified data points, if available.

% Copyright (c) 2002-2005 Finn Lindgren
% $Id: normmix_classify.m 2104 2005-10-31 15:09:32Z finn $

% Parse input parameters
if nargin>1
  if (nargin<4), gt = []; end
  p = normmix_posterior(x,theta,prior,gt);
else
  p = x;
end

% Classify:
[tmp,cl] = max(p,[],2);

[n,K] = size(p);
ind = 1:K;
cl_ind = (ind(ones(n,1),:)==cl(:,ones(1,K)))*1; % *1 to avoid logical array.
