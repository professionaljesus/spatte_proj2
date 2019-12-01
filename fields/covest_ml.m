function [par,beta] = covest_ml(D, z, covf, par0, X, method, par_start)
% covest_ml estimates parameters of a covariance using Maximum Likelihood.
%
% par = covest_ml(D, z)
%  or
% par = covest_ml(D, z, covf, par_fixed)
%  or
% [par,beta]=covest_ml(D, y, covf, par_fixed,X)
%  or
% [par,beta]=covest_ml(D, y, covf, [], X)
%  or
% [par,beta]=covest_ml(D, y, covf, par_fixed, X, method, par_start)
%
% D - The distance matrix for the observations, or coordinates. If D is not
%     a square matrix then covest_ml computed a distance matrix by calling:
%     D=distance_matrix(D).
% z - are the residuals  z=y-mu,  for known expectation  mu.
% y - are the observations, if mu is unknown, but X supplied.
%
% Optional parameters:
% covf - Covariance function to use, one of the following strings:
%        'matern' (default), 'cauchy', 'exponential', 'gaussian',
%        'spherical'
% par_fixed - Has the same structure as  par.  (default=zeros(n,1))
%             Allows setting some parameters to known values.
%             A non-zero entry for a parameter fixes that parameter
%             to that value.
% X - The regression/expectation basis matrix X. If X is supplied, the
%     optimisation also estimates beta, for mu=X*beta. If X is supplied,
%     beta are the estimated expectation parameters. 
% method - Should we use standard ML (profile likelihood) or REML. Only
%          matters if X is given. one of the following strings: 
%          'ml' (default), 'reml'
% par_start - A vector of initial values for the optimisation.
%
% Returns  estimates of 
%   par = [par_cov; sigma2_epsilon] - Order of parameters in par_cov is the
%          same as order in the corresponding covariance function, e.g. 
%          par_cov = [sigma2;kappa;nu] for matern_covariance.
%   beta - regression parameters for the case of an unknown mean.


% $Id: covest_ml.m 5090 2017-11-05 19:35:05Z johanl $

if nargin<3 || isempty(covf), covf='matern'; end
if nargin<4, par0 = []; end
if nargin<5, X = []; end
if nargin<6 || isempty(method), method='ml'; end
if nargin<7, par_start = []; end

%Check which covariance function we want
covf = validatestring(covf, {'matern', 'cauchy', 'exponential', ...
	'gaussian', 'spherical'}, 'covest_ml', 'covf', 3);

%Check which method we want (ml for the X=[] case)
method = validatestring(method, {'ml', 'reml'}, 'covest_ml', 'method', 6);
if isempty(X), method='ml'; end
REML = strcmp(method, 'reml');

%do we need to compute the distance matrix
if size(D,1)~=size(D,2), D = distance_matrix(D); end
%find unique distances
[dunique,I,J] = unique(D(D>0));

%rough estimate fo field variance
if isempty(X), s2 = var(z); else, s2 = var(z-X*(X\z)); end

[covf, par0, par_start] = covest_init(covf, par0, par_start, max(dunique), s2);

%indicator of fixed parameters
fixed = par0~=0;

%define loss function that accounts for fixed parameters
neg_log_like = @(par) covest_ml_loss(D, z, exp(par).*(~fixed) + par0, ...
	X, I, J, covf, REML) + (fixed'*(exp(par)-par0).^2);
%minimize loss function (negative log-likelihood)
par = fminsearch(neg_log_like, log(par_start));
%extract estimated parameters
par = exp(par).*(~fixed) + par0;

if ~isempty(X) % X supplied, estimate beta
	[~,beta] = covest_ml_loss(D, z, par, X, I, J, covf, REML);
else
	beta=[];
end

%% loss function for the covariance function
function [nlogf, beta] = covest_ml_loss(D, z, par, X, I, J, covf, REML)
%extract parameters
par_covf = par(1:end-1);
sigma2_eps = par(end);
%compute covariance function
Sigma_yy = covf(D,par_covf,I,J) + sigma2_eps*eye(length(z));
%compute Choleskey factor
[R,p] = chol(Sigma_yy);
%set default beta value
beta = [];
%Choleskey fail -> return large value
if p~=0
	nlogf = realmax/2;
	return; 
end

if isempty(X)
  nlogf = sum(log(diag(R))) + 0.5*sum((z'/R).^2);
else
  tmp = ([z,X]'/R)';
  beta = tmp(:,2:end)\tmp(:,1);
  nlogf = sum(log(diag(R))) + 0.5*sum((tmp(:,1)-tmp(:,2:end)*beta).^2);
	if REML
		XiSX = tmp(:,2:end)'*tmp(:,2:end);
		%compute Choleskey factor
		[R,p] = chol(XiSX);
		%Choleskey fail -> return large value
		if p~=0
			nlogf = realmax/2;
			return;
		end
		nlogf = nlogf + sum(log(diag(R)));
	end
end
