function theta_hat=gmrf_param_map(gmrf_negloglike,theta_start)
% GMRF_PARAM_MAP  Finds the MAP estimate in a simple field model.
%
%  theta_hat=gmrf_param_map(gmrf_negloglike,theta_start)
%
% theta = [log([tau2; kappa; sigma2_epsilon]); beta]
% gmrf_negloglike = a function handle to a function calculating the
%                   negated log-likelihood, see gmrf_negloglike_skeleton
%
% Example, with data-adapted starting values:
%  kappa_0 = sqrt(8)/1;
%  R_tmp = chol(kappa_0^4*C+2*kappa_0^2*G+G2);
%  tau2_0 = var(Y)*R_tmp(end,end)^2;
%  beta_0 = (A*B)\Y;
%  theta_start = [log([tau2_0;kappa_0;0.01*tau2_0]);beta_0];
%  theta_hat = gmrf_param_map(@(th) gmrf_negloglike(th,Y,A,C,G,G2,B),...
%                             theta_start)
%
% Note: Currently, this function simply calls fminunc; Feel free to
%       call fminunc or fminsearch by yourself if you want to!

% $Id: gmrf_param_map.m 4836 2014-12-10 11:09:32Z johanl $

theta_hat = fminunc(@(th) gmrf_negloglike(th),theta_start);
