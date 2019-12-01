function output=gmrf_mcmc_skeleton(gmrf_negloglike,theta_start,Sigma_propose,T,extra)
% GMRF_MCMC  Simulate GMRF parameters with MCMC, and perform calculations
%
%  output = gmrf_mcmc(theta_hat,gmrf_negloglike,extra)
%
% gmrf_negloglike = a function handle to a function calculating the
%                   negated log-likelihood, see gmrf_negloglike_skeleton
% theta_start = the starting point for the simulation,
%               [log([tau2; kappa; sigma2_epsilon]); beta]
% Sigma_propose = the proposal covariance, see lecture F07, 2009.
% T = the number of samples to generate
% extra = any extra variables you need for your calculations
%
% Skeleton example:
%  theta_hat = ...;
%  Sigma_propose = 2.4/6*inv(gmrf_param_hessian(...))
%  gmrf_mcmc_skeleton(@(th) gmrf_negloglike(th,y,A,C,G,G2,B),...
%                     theta_hat,Sigma_propose,1000,...)
%
% (Note: Since the simulation should be started in the MAP estimate of
%        theta, we ignore the otherwise needed burn-in period.)
%
% This is only a skeleton for Home Assignment 2.

% $Id: gmrf_mcmc_skeleton.m 5107 2017-11-12 13:35:17Z johanl $

% Remove this line from your copy:
warning(['This is only a skeleton function!  Copy it and fill in the blanks!'])

theta_t = theta_start;
negloglike_t = gmrf_negloglike(theta_t);

R_propose = chol(Sigma_propose);

for t=1:T
  theta_prop = theta_t + R_propose'*randn(length(theta_start),1);
  negloglike_prop = gmrf_negloglike(theta_prop);
  if (log(rand(1,1)) <= -negloglike_prop+negloglike_t)
    % Accept the proposal:
    theta_t = theta_prop;
    negloglike_t = negloglike_prop;
  end

  tau2 = exp(theta_t(1));
  kappa = exp(theta_t(2));
  sigma2eps = exp(theta_t(3));
  beta = theta_t(4:end);

  output.theta(:,t) = theta_t; % Save the theta-values, if you want to.

  % Fill in your own calculations for the current theta_t here:
  something(t) = 0;
end

% Store the result in the output variable:
output.something = 0;
