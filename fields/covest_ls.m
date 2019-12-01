function par=covest_ls(rhat, s2hat, m, n, d, covf, par0, par_start)
% covest_ls estimates parameters of a covariance using Least Squares.
%
% par=covest_ls(rhat, s2hat, m, n, d, covf, par_fixed)
%
% par = [sigma2; kappa; nu; sigma2_epsilon]
% rhat,s2hat,m,n,d  are the outputs from covest_nonparametric
%
% Optional parameters:
% covf - Covariance function to use, one of the following strings:
%        'matern' (default), 'cauchy', 'exponential', 'gaussian',
%        'spherical'
% par_fixed - Has the same structure as  par.  (default=zeros(n,1))
%             Allows setting some parameters to known values.
%             A non-zero entry for a parameter fixes that parameter
%             to that value.
% par_start - A vector of initial values for the optimisation.
%
% Returns  estimates of 
%   par = [par_cov; sigma2_epsilon] - Order of parameters in par_cov is the
%          same as order in the corresponding covariance function, e.g. 
%          par_cov = [sigma2;kappa;nu] for matern_covariance.

% $Id: covest_ls.m 5090 2017-11-05 19:35:05Z johanl $

if nargin<6 || isempty(covf), covf='matern'; end
if nargin<7, par0 = []; end
if nargin<8, par_start = []; end

%Check which covariance function we want
covf = validatestring(covf, {'matern', 'cauchy', 'exponential', ...
	'gaussian', 'spherical'}, 'covest_ls', 'covf', 6);

[covf, par0, par_start] = covest_init(covf, par0, par_start, max(d), s2hat);

%indicator of fixed parameters
fixed = par0~=0;

%define loss function that accounts for fixed parameters
WLS_loss = @(par) covest_ls_loss(rhat, s2hat, m, n, d, ...
	exp(par).*(~fixed) + par0, covf) + (fixed'*(exp(par)-par0).^2);
%minimize loss function (negative log-likelihood)
par = fminsearch(WLS_loss, log(par_start));
%extract estimated parameters
par = exp(par).*(~fixed) + par0;


function S=covest_ls_loss(rhat, s2hat, m, n, d, par, covf)
%extract parameters
par_covf = par(1:end-1);
sigma2_eps = par(end);
%check for non empty bins
ok = ~isnan(rhat);
%compute covariance function
r = covf(d(ok),par_covf,[],[]);
%compute WLS
S = n*(s2hat - (par_covf(1)+sigma2_eps))^2 + ...
	sum( m(ok).*(rhat(ok)-r).^2 );
