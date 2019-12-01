function [theta0,prior0]=normmix_kmeans(x,K,maxiter,plotflag)
% NORMMIX_KMEANS Use a single K-means to get a random inital estimate
%                of the parameters in a Gaussian mixture model.
%
% [theta0,prior0]=normmix_kmeans(x,K)
% [theta0,prior0]=normmix_kmeans(x,K,maxiter,plotflag)
%
% x: n-by-d matrix
% theta0{k}.mu: 1-by-d matrix, class expected value.
% theta0{k}.Sigma: d-by-d matrix, class covariance.
% prior0: 1-by-K matrix, the class probabilities
%
% maxiter: the maximum number of iterations in K-means (default=1)
%
% plotflag==0: no plots. (default)
% plotflag>=1: plots the succesive resulting classifications.
% plotflag>=2: plots the parameter estimate trails.
% plotflag>=3: plots convergence monitoring information.

% Copyright (c) 2002-2005 Finn Lindgren
% $Id: normmix_kmeans.m 4586 2012-10-08 16:18:33Z johanl $

% Parse input parameters
if nargin<3, maxiter = []; end
if nargin<4, plotflag = []; end
if isempty(maxiter), maxiter = 1; end
if isempty(plotflag), plotflag = 0; end

%% plot %%
if plotflag
  figure(1),clf
end
%% plot %%

% Use some steps of K-means for rough initial estimates:
[cl,theta0] = kmeans(x,K,maxiter,plotflag);
% pi-estimates (use uniform pi):
prior0 = ones(1,K)/K;
% Sigma-estimates (use common, isotropic Sigma):
n = size(x,1);
d = size(x,2);
Sigma = 0;
for k=1:K
  y = x(cl==k,:)-repmat(theta0{k}.mu,[sum(cl==k),1]);
  Sigma = Sigma + sum(sum(y.*y,1),2);
end
Sigma = eye(d)*Sigma/n/d;
for k=1:K
  theta0{k}.Sigma = Sigma;
end
