function [theta,prior,p,samples]=normmix_gibbs(x,K,Nsim,plotflag,prior0)
% NORMMIX_GIBBS Sample parameters in a Gaussian mixture model using Gibbs.
%
% [theta,prior]=normmix_gibbs(x,K)
% [theta,prior]=normmix_gibbs(x,K,Nsim)
% [theta,prior]=normmix_gibbs(x,K,Nsim,plotflag)
% [...]=normmix_gibbs(x,theta0,...)
% [...]=normmix_gibbs(x,theta0,...,prior0)
% [...]=normmix_gibbs(x,cl0,...)
% [...,p]=normmix_gibbs(...)
% [...,p,samples]=normmix_gibbs(...)
%
% x: n-by-d matrix
%
% The following will be posterior means from the Gibbs sampling.
% theta{k}.mu: 1-by-d matrix, class expected value.
% theta{k}.Sigma: d-by-d matrix, class covariance.
% prior: 1-by-K matrix, the class probabilities
% p: The posterior class probabilities, n-by-K matrix.
%
% Nsim: Nsim(1) Total number of Gibbs-samples
%       Nsim(2) Number of Gibbs-samples used for burn-in
%               Default: Nsim = [200, 100]
% samples: A structure with Gibbs samples
%    samples.mu: Nsim(1)-by-d-by-K matrix
%    samples.p: Nsim(1)-by-K matrix
%
% plotflag==0: no plots. (default)
% plotflag>=1: plots the succesive resulting classifications.
%              For d==1, plots empirical and model densities,
%                        and the posterior probabilities.
% plotflag>=2 also plot densities for the Gaussian mixture model
%             that would generate the same classifications,
%             and the posterior class probabilities.
%             (only for 1D-data)
% plotflag>=3: plots the parameter estimate trails.
%
% theta0: the initial class parameter estimates.
% prior0: 1-by-K matrix, the initial class probability estimates.
% cl0: vector of initial classifications used to estimate intial
%      parameters.

% $Id: normmix_gibbs.m 5038 2016-12-05 00:14:04Z johanl $

% Parse input parameters
if nargin<3, Nsim = []; end
if nargin<4, plotflag = []; end
if nargin<5, prior0 = []; end
if isempty(Nsim), Nsim = [200 100]; end
if (length(Nsim)<2), Nsim(2) = round(Nsim/2); end
if Nsim(1)<=Nsim(2)
  error('Nsim(1)<=Nsim(2) - Only Burn-in specified'); 
end
if isempty(plotflag), plotflag = 0; end

%% plot %%
if plotflag
  if plotflag>=3, figure(2),clf, end
  figure(1),clf
end
%% plot %%

%input sizes
[n,d] = size(x);
if iscell(K)  % Initial estimates were supplied.
  theta0 = K;
  K = length(theta0);
  if isempty(prior0), prior0 = ones(1,K)/K; end
elseif isvector(K) && size(K,1)==size(x,1)
  cl = K;
  K = max(cl);
  theta0 = cell(1,K);
  for i=1:K
    theta0{i}.mu = mean(x(cl==i,:),1);
    theta0{i}.Sigma = diag(var(x(cl==i,:),[],1));
  end
  prior0 = ones(1,K)/K;
elseif isscalar(K)
  % No initial theta estimates available.
  if isempty(prior0)
    [theta0,prior0] = normmix_kmeans(x,K,1,plotflag);
  else
    theta0 = normmix_kmeans(x,K,1,plotflag);
  end
else
  error('Unknown input type for K')
end

%storage of samples (past burnin)
p = zeros(n,K);
prior = zeros(1,K);
theta = cell(1,K);
for k=1:K
  theta{k}.mu = zeros(1,d);
  theta{k}.Sigma = zeros(d,d);
end
Nsample = 0;
n_k = zeros(1,K);

%should we save the sample-paths
if nargout>3 || plotflag
  samples.mu = zeros(Nsim(1), d, K);
  samples.p = zeros(Nsim(1), K);
end

for i=1:Nsim(1)
  %% Compute posterior probabilities %%
  p_tmp = normmix_posterior(x, theta0, prior0);
  %% Sample class belongings from the posterior %%
  z = sum(bsxfun(@gt, rand(size(p_tmp,1),1), cumsum(p_tmp,2)),2)+1;
  %% Sample new class parameters %%
  for k=1:K
    Ind = (z==k);
    [theta0{k}.mu, theta0{k}.Sigma] = gibbs_mu_sigma(x(Ind,:));
    n_k(k) = sum(Ind);
  end
  %% Sample new class probabilities %%
  prior0 = dirichletrnd(n_k+1);
  
  %% if we're past burn-in save samples %%
  if i>Nsim(2)
    p = p + p_tmp;
    prior = prior + prior0;
    for k=1:K
      theta{k}.mu = theta{k}.mu + theta0{k}.mu;
      theta{k}.Sigma = theta{k}.Sigma + theta0{k}.Sigma;
    end
    Nsample = Nsample+1;
  end
  
  %% save sample paths %%
  if nargout>3 || plotflag
    for k=1:K
      samples.mu(i,:,k) = theta0{k}.mu;
    end
    samples.p(i,:) = prior0;
  end

  %% plot %%
  if plotflag
    doplot(x, z, theta0, prior0, p_tmp, plotflag);
    if size(x,2)>=3
      hold on
      plot3(squeeze(samples.mu(1:i,1,:)),...
	    squeeze(samples.mu(1:i,2,:)),...
	    squeeze(samples.mu(1:i,3,:)),'k-')
      hold off
    elseif size(x,2)>=2
      hold on
      plot(squeeze(samples.mu(1:i,1,:)),squeeze(samples.mu(1:i,2,:)),'k-')
      hold off
    end
    title(i)
    if plotflag>=3
      figure(2)
      subplot(211)
      plot(1:i,samples.p(1:i,:))
      title('\pi')
      subplot(212)
      plot(1:i,squeeze(samples.mu(1:i,1,:)))
      title('\mu_1')
      figure(1)
    end
    drawnow
  end
end

%% Compute posterior means %%
p = p / Nsample;
prior = prior / Nsample;
for k=1:K
  theta{k}.mu = theta{k}.mu / Nsample;
  theta{k}.Sigma = theta{k}.Sigma / Nsample;
end

%%%%%
function doplot(x,cl,theta,prior,p,plotflag)
col = 'bgrcmyk';
n = size(x,1);
d = size(x,2);
K = length(theta);
if (d==1)
  bins = linspace(min(x),max(x),ceil(sqrt(n)));
  if (plotflag>=2)
    subplot(211)
  else
    subplot(111)
  end
  N = cell(1,K);
  p_model = zeros(K, length(bins));
  for k=1:K
    N{k} = histc(x(cl==k),bins)';
    stairs(bins,N{k}./gradient(bins)/n,col(k))
    hold on
    y = bins-theta{k}.mu;
    p_model(k,:) = prior(k)*exp(-y.^2/(2*theta{k}.Sigma)) / ...
        ((2*pi*theta{k}.Sigma)^(1/2));
    plot(bins,p_model(k,:),col(k))
  end
  plot(bins,sum(p_model,1),'k')
  hold off
  if (plotflag>=2)
    subplot(212)
    for k=1:K
      plot(x,p(:,k),['.' col(k)])
      hold on
    end
    hold off
  end
elseif (d>=3)
  for k=1:K
    plot3(x(cl==k,1),x(cl==k,2),x(cl==k,3),['.' col(k)],...
	 theta{k}.mu(1),theta{k}.mu(2),theta{k}.mu(3),['o' col(k)]);
    hold on
  end
  hold off
else % d==2 
  for k=1:K
    plot(x(cl==k,1),x(cl==k,2),['.' col(k)],...
	 theta{k}.mu(1),theta{k}.mu(2),['o' col(k)]);
    hold on
  end
  hold off
end
