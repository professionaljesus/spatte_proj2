function negloglike=gmrf_negloglike_skeleton(theta,y,A,C,G,G2,B)
% GMRF_NEGLOGLIKE_SKELETON  Calculate the GMRF data likelihood, x integrated out
%
% negloglike=gmrf_negloglike_skeleton(theta,y,A,C,G,G2,B)
%
% theta = [log([tau2; kappa; sigma2_epsilon]); beta]
% y = the data vector, as a column with n elements
% A = the observation matrix, sparse n-by-N
% C,G,G2 = matrices used to build a Matérn-like precision,
%          see matern_prec_matrices, sparse N-by-N
% B = the expectation basis matrix, N-by-length(beta)
%
% This is only a skeleton for Home Assignment 2.

% $Id: gmrf_negloglike_skeleton.m 5107 2017-11-12 13:35:17Z johanl $

n = size(y,1);
tau2 = exp(theta(1));
kappa = exp(theta(2));
sigma2eps = exp(theta(3));
beta = theta(4:end);

% Remove this line from your copy:
warning('This is only a skeleton function!  Copy it and fill in the blanks!')

mu = [];
Q_x = [];

%compute optimal reordering
order = amd(Q_x);
%reorder relevant matrices.
Q_x = Q_x(order,order);
mu = mu(order);
A = A(:,order);

%Compute choleskey factor
[R_x,p_x] = chol(Q_x);

if p_x~=0
  %choleskey factor fail -> (almost) semidefinite matrix -> 
  %-> det(Q) ~ 0 -> log(det(Q)) ~ -inf -> negloglike ~ inf
  %Set negloglike to a REALLY big value
  negloglike = realmax;
else
  Q_xy = [];
  R_xy = chol(Q_xy);
  mu_xy = [];

  negloglike = 0;
end

%if we want to use mu_xy we need to invert the reordering
%mu_xy(order) = mu_xy;
