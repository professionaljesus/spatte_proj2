function [negloglike,mu_xy,Q_xy]=gmrf_negloglike(theta, y, A, AtA, spde, N_beta)
% GMRF_NEGLOGLIKE  Calculate the GMRF data likelihood
%
% negloglike=gmrf_negloglike_skeleton(theta, y, A, AtA, spde)
%
% theta = [tau2_1; tau2_2; tau2_3; kappa_1; kappa_2; kappa_3; 
%               sigma2_epsilon]);]
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-(3*N)
% AtA = precomputed A'*A, sparse 3*N-by-3*N
% spde = struct with the three matrices (C,G,G2) used to build a Matern-like
%        precision, sparse N_x-by-N_x. (and N=N_x+N_beta)
% N_beta = vector with the number of covariates for the expectation in each
%          of the the three fields
%
% This is only a skeleton for Home Assignment 3.

% $Id: gmrf_negloglike.m 4424 2011-09-19 14:58:51Z johanl $

%number of observations
n = size(y,1);

%extract parameters for the three Q matrices
tau = exp( theta(1:3) );
kappa2 = exp( theta(4:6) );
%and the nugget
sigma2eps = exp( theta(7) );

%create Qbeta (same prior for all betas, but possibly different number of
%covariates, i.e. different sizes)
Qbeta = cell(3,1);
for i=1:3
  Qbeta{i} = 1e-6 * speye(N_beta(i));
end

%create the three Q matrices (one for each time-trend)
Q = cell(3,1);
for i=1:length(Q)
  Q{i} = []; %using tau(i) and kappa2(i)
end

%combine the Q matrices to one large matrix
Q_all = blkdiag(Q{1}, ...?

%and construct the posterior precision
Q_xy = [];

%compute Choleskey factor for each of the Q{i} matrices
p_chol = zeros(3,1);
R = cell(length(Q),1);
for i=1:length(R)
  [R{i}, p_chol(i)] = chol(Q{i});
end
%and compute choleskey for posterior matrix (here the reordering matters,
%fell free to examin performance WITHOUT the reordering!)
p = amd(Q_xy);
[R_xy,pxy] = chol( Q_xy(p,p) );
if any(p_chol~=0) || pxy~=0
  %choleskey factor fail -> (almost) semidefinite matrix -> 
  %-> det(Q) ~ 0 -> log(det(Q)) ~ -inf -> negloglike ~ inf
  %Set negloglike to a REALLY big value
  negloglike = realmax;
  %and set posterio mean as "bad" (matlab needs a return value)
  mu_xy = nan;
  return;
end

%compute posterior mean
%first compute b
b_xy = [];
%reorder
b_xy = b_xy(p);
%compute remainder of the mean (use R_xy)
mu_xy = [];
%and revert the order
mu_xy(p) = mu_xy;

%compute negative loglikelihood (first half of Lecture 7)
negloglike = [];

%output for analysis of optimisation run
disp( reshape(theta(:), [], numel(theta)) )
disp(negloglike)