function [theta,prior,p,converged]=normmix_em(x,K,convergence,plotflag,prior0,gt)
% NORMMIX_EM Estimate parameters in a Gaussian mixture model.
%
% [theta,prior]=normmix_em(x,K)
% [theta,prior]=normmix_em(x,K,convergence)
% [theta,prior]=normmix_em(x,K,convergence,plotflag)
% [...]=normmix_em(x,theta0,...)
% [...]=normmix_em(x,theta0,...,prior0)
% [...,p]=normmix_em(...)
% [...,p,converged]=normmix_em(...)
% [...]=normmix_em(...,gt) % Ground truth data, relies on support in
%                            normmix_posterior
%
% x: n-by-d matrix
% theta{k}.mu: 1-by-d matrix, class expected value.
% theta{k}.Sigma: d-by-d matrix, class covariance.
% prior: 1-by-K matrix, the class probabilities
% p: The posterior class probabilities, n-by-K matrix.
%
% convergence: convergence(1) is the maximum number of iterations
%              convergence(2) is the p-difference tolerance
%              Default: convergence = [200,5e-4]
% converged: 1 if the maximal w-difference was < convergence(2) in
%            the final iteration, 0 otherwise.
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
% plotflag>=4: plots convergence monitoring information.
%
% theta0: the initial class parameter estimates.
% prior0: 1-by-K matrix, the initial class probability estimates.
% gt:     ground truth, see NORMMIX_POSTERIOR

% Copyright (c) 2002-2005 Finn Lindgren
% $Id: normmix_em.m 4590 2012-10-08 20:37:55Z johanl $

% Parse input parameters
if nargin<3, convergence = []; end
if nargin<4, plotflag = []; end
if nargin<5, prior0 = []; end
if nargin<6, gt = []; end
if isempty(convergence), convergence = 200; end
if (length(convergence)<2), convergence(2) = 5e-4; end
if isempty(plotflag), plotflag = 0; end

%% plot %%
if plotflag
  if plotflag>=4, figure(3),clf, end
  if plotflag>=3, figure(2),clf, end
  figure(1),clf
end
%% plot %%

if iscell(K)  % Initial estimates were supplied.
  theta0 = K;
  K = length(theta0);
  if isempty(prior0), prior0 = ones(1,K)/K; end
else % No initial theta estimates available.
  if ~isempty(gt)
    warning('fms150:gtHasNoEffectInNormmixKmeans',...
            ['Ground truth available, but normmix_kmeans ',...
             'can not not use it!\n',...
             'Compute initial parameter estimates yourself!'])
  end
  if isempty(prior0)
    [theta0,prior0] = normmix_kmeans(x,K,1,plotflag);
  else
    theta0 = normmix_kmeans(x,K,1,plotflag);
  end
end

[n,d] = size(x);
prior = prior0;
theta = theta0;

% Perform the E-step:
p = E_step(x,theta,prior,gt);

% Classify, for plotting:
[tmp,cl] = max(p,[],2);

%% plot %%
if plotflag
  doplot(x,cl,theta,prior,p,plotflag);
end
for k=1:K
  prior_history(1,k,1) = prior(k);
  mu_history(k,:,1) = theta{k}.mu;
  Sigma_history(:,:,k,1) = theta{k}.Sigma;
end
%% plot %%

loop = 0;
done = 0;
while (~done)
  loop = loop+1;
  % For the stopping criterion:
  p_old = p;
  % The E-step:
  if (loop>1) % Work not already done above.
    p = E_step(x,theta,prior,gt);
    
    % Classify, for plotting:
    [tmp,cl] = max(p,[],2);
  end

  % The M-step:
  % New pi-estimates:
  prior = sum(p,1)/n;

  % New mu and Sigma estimates:
  for k=1:K
    theta{k}.mu = sum(x.*repmat(p(:,k),[1,d]))/sum(p(:,k));
    y = x-repmat(theta{k}.mu,[n,1]);
    theta{k}.Sigma = y'*(y.*repmat(p(:,k),[1,d]))/sum(p(:,k));
  end

  % The stopping criterion:
  p_diff(loop) = max(abs((p_old(:)-p(:))));
  converged = (loop>1) & (p_diff(loop)<convergence(2));
  done = converged | (loop>=convergence(1));

  %% plot %%
  if plotflag
    for k=1:K
      prior_history(1,k,loop+1) = prior(k);
      mu_history(k,:,loop+1) = theta{k}.mu;
      Sigma_history(:,:,k,loop+1) = theta{k}.Sigma;
    end
    doplot(x,cl,theta,prior,p,plotflag);
    if size(x,2)>=3
      hold on
      plot3(squeeze(mu_history(:,1,:))',...
	    squeeze(mu_history(:,2,:))',...
	    squeeze(mu_history(:,3,:))','k-')
      hold off
    elseif size(x,2)>=2
      hold on
      plot(squeeze(mu_history(:,1,:))',squeeze(mu_history(:,2,:))','k-')
      hold off
    end
    title(loop)
    if (plotflag>=3)
      figure(2)
      subplot(311)
      plot(0:loop,squeeze(prior_history)')
      title('\pi')
      subplot(312)
      plot(0:loop,squeeze(mu_history(:,1,:))')
      title('\mu_1')
      subplot(313)
      plot(0:loop,squeeze(Sigma_history(1,1,:,:))')
      title('\Sigma_{11}')
      if (plotflag>=4)
        if (loop>1)
          figure(3)
          subplot(211)
          semilogy(1:loop,p_diff)
          title('Maximal p-difference')
          subplot(212)
          lim = p_diff(loop);
          p_diff_n = histc(p_old(:)-p(:),...
                           linspace(-lim,lim,100))*100/length(p(:));
          bar(linspace(-lim,lim,100),p_diff_n.^0.25,'histc');
          axis([-lim,lim,0,max(p_diff_n).^0.25])
          title('p-difference histogram')
        end
      end
      figure(1)
    end
    drawnow
  end
  %% plot %%
end

%%%%%
function p=E_step(x,theta,prior,gt)

% Calculate the posterior class probabilities p:
p = normmix_posterior(x,theta,prior,gt);


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
