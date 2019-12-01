function [cl,theta]=kmeans(x,K,maxiter,plotflag)
% KMEANS Classify data using the K-means algorithm.
%
% [cl,theta]=kmeans(x,K)
% [cl,theta]=kmeans(x,K,maxiter)
% [cl,theta]=kmeans(x,K,maxiter,plotflag)
% [cl,theta]=kmeans(x,theta0,...)
% [cl,theta]=kmeans(x,cl0,...)
%
% K: The number of clusters
% theta0: Initial cluster centers as theta0{k}.mu
% cl0: Initial classification, used to compute initial cluster centers as
%      theta0{k}.mu = mean(x(cl==k,:),1);
% maxiter: The maximum number of iterations. (default=Inf)
% plotflag: ==0 for no plotting (default)
%           >=1 for d>=2, plot scatter plot of data, with
%               colour indicating classifications
%               for d==1, plot empirical densities
%           >=2 also plot densities for the Gaussian mixture model
%               that would generate the same classifications,
%               and the posterior class probabilities.
%               (only for 1D-data)
%
% x: n-by-d matrix
% theta{k}.mu: 1-by-d matrix

% Copyright (c) 2002,2003 Finn Lindgren
% $Id: kmeans.m 5038 2016-12-05 00:14:04Z johanl $

if nargin<3, maxiter = []; end
if nargin<4, plotflag = []; end
if isempty(maxiter), maxiter = inf; end
if isempty(plotflag), plotflag = 0; end

[n,d] = size(x);

if iscell(K)
  theta = K;
  K = numel(K);
  mu = zeros(K,d);
  for i=1:K, mu(i,:) = theta{K}.mu; end
  %check that mu's are unique
  if size(unique(mu,'rows'),1)<K
    error('class starting points are not unique')
  end
elseif isvector(K) && size(K,1)==size(x,1)
  cl = K;
  K = max(cl);
  mu = nan(K,d);
  for i=1:K, mu(i,:)=mean(x(cl==i,:),1); end
  if any(isnan(mu))
    error('not all classes in 1..K included in K=cl vector')
  end
elseif isscalar(K)
  % Find unique starting points:
  if (n<K), error('Not enough data!'), end
  start_idx = ceil(rand(K,1)*n);
  mu = x(start_idx,:);
  while (size(unique(mu,'rows'),1)<K)
    start_idx = ceil(rand(K,1)*n);
    mu = x(start_idx,:);
  end
else
  error('Unknown input type for K')
end


cl_old = ones(n,1);
cl = zeros(n,1);

%% plot %%
if plotflag
  doplot(x,cl,mu,plotflag); drawnow
end
%% plot %%

loop = 0;
mu_history(:,:,loop+1) = mu;
dist2 = zeros(size(x,1),K);
while (~all(cl==cl_old)) && (loop<maxiter)
  loop = loop+1;
  cl_old = cl;
  % Squared distances:
  for i=1:K
    dist2(:,i) = sum(bsxfun(@minus, x, mu(i,:)).^2,2);
  end
  % Classify:
  [~,cl] = min(dist2,[],2);
  % New mu estimates
  for i=1:K
    mu(i,:) = mean(x(cl==i,:),1);
  end

  %% plot %%
  if plotflag
    mu_history(:,:,loop+1) = mu;
    doplot(x,cl,mu,plotflag);
    if size(x,2)>1
      hold on
      plot(squeeze(mu_history(:,1,:))',squeeze(mu_history(:,2,:))','k-')
      hold off
    end
    title(loop)
    drawnow
  end
  %% plot %%
end

% Collect output:
theta = cell(K,1);
for k=1:K
  theta{k}.mu = mu(k,:);
end

%%%%%
function doplot(x,cl,mu,plotflag)
col = 'bgrcmyk';
n = size(x,1);
d = size(x,2);
K = size(mu,1);
if (d==1)
  if (any(cl>0))
    if (plotflag>=2)
      subplot(211)
    else
      subplot(111)
    end
    bins = linspace(min(x),max(x),ceil(sqrt(n)));
    N = cell(K,1);
    for k=1:K
      N{k} = histc(x(cl==k),bins)';
      stairs(bins,N{k}./gradient(bins)/n,col(k))
      hold on
    end
    if (plotflag>=2)
      prior = ones(1,K)/K;
      Sigma = 0;
      for k=1:K
        y = x(cl==k,:)-mu(k);
        Sigma = Sigma + y'*y;
      end
      Sigma = Sigma/size(x,1);
      p = zeros(K,length(bins));
      for k=1:K
        y = bins-mu(k);
        p(k,:) = prior(k)*exp(-y.^2/(2*Sigma)) / ...
                 ((2*pi*Sigma)^(1/2));
        plot(bins,p(k,:),col(k))
      end
      plot(bins,sum(p,1),'k')
      hold off
      subplot(212)
      for k=1:K
        plot(bins,p(k,:)./sum(p,1),col(k))
        hold on
      end
      hold off
    end
  end
else
  min1 = min(x(:,1));
  max1 = max(x(:,1));
  min2 = min(x(:,2));
  max2 = max(x(:,2));
  w = max1-min1;
  h = max2-min2;
  c = cos(linspace(0,2*pi,10))*w*100+(min1+max1)/2;
  s = sin(linspace(0,2*pi,10))*h*100+(min2+max2)/2;
  voronoi([mu(:,1);c'],[mu(:,2);s'])
  hold on
  for k=1:K
    plot(x(cl==k,1),x(cl==k,2),['.' col(k)],...
	 mu(k,1),mu(k,2),['o' col(k)]);
  end
  hold off
  axis([[min1,max1]+[-1,1]*w*0.1,[min2,max2]+[-1,1]*h*0.1])
end
