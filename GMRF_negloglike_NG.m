function negloglike = GMRF_negloglike_NG(theta, model_order, y, A, B, G, E)
% GMRF_NEGLOGLIKE_NG_SKELETON Calculate the GMRF likelihood for non-Gaussian observations
%
% negloglike = GMRF_negloglike_NG(theta, y, A, B, G, E, model_order)
%
% theta - log of parameters
% y - data vector, as a column with n elements
% A - Observation matrix, sparse n-by-N
% B - Matrix of covariates, matrix of size n-by-Nbeta 
% C,G,G2 = matrices used to build a Matern-like CAR/SAR precision,
%          see matern_prec_matrices, sparse N-by-N
% qbeta = Precision for the regression parameters, scalar.
% E = The population count in each region, used for Poisson observations
%
% This is only a skeleton for Home Assignment 2.

% $Id: gmrf_negloglike_NG_skeleton.m 5109 2017-11-12 20:08:26Z johanl $

% Remove this line from your copy:
% error('This is only a skeleton function!  Copy, rename, and fill in the blanks!')

%ensure that E=1 if E not given (i.e. same/no population weight in all regions)
if nargin<7, E=1; end

%extract parameters (and transform to natural parameters)
tau = exp(theta(1));

%compute Q matrices  (for an intrinsic CAR(1) or SAR(1) process)
% help matern_prec_matrices: Q_x = kappa^4*C + 2*kappa^2*G + G2 ... but kappa = 0 in our case.
if model_order == 2
    Q_x = tau.*G*G; % SAR(1)
elseif model_order == 1
    Q_x = tau.*G; % CAR(1)
else
    disp('Unvalid model order. Model order must be 1 or 2')
end

%compute Q for beta-parameters
q_beta = 1e-3;
Qbeta = q_beta.*speye(size(B, 2));

% Create different Aall and Qall depending on number of parameters fed in 
% to the function. 
% Also compute the determinant for Qall here. Reordering of Qall does not 
% matter for the determinant in our case.
if length(theta) == 2
    %compute Q for nugget-parameters
    q_eps = exp(theta(2));
    Qeps = q_eps.*speye(size(A, 2));
    
    %combine all components of Q using blkdiag
    Qall = blkdiag(Q_x, Qeps, Qbeta);
    
    % also compute determinante of Qall (might be simplified)
    Qall_det = log(tau)*length(Q_x) + log(q_beta)*length(Qbeta) + log(q_eps)*length(Qeps);
    
    %also compute the observation matric by combining A and B matrices
    Aall = [A A B];
    
elseif length(theta) == 1
    %combine all components of Q using blkdiag
    Qall = blkdiag(Q_x, Qbeta);
    
    % also compute determinante of Qall (might be simplified)
    Qall_det = log(tau)*length(Q_x) + log(q_beta)*length(Qbeta);
    
    %also compute the observation matric by combining A and B matrices
    Aall = [A B];
else
    disp('Unvalid dimension of parameter vector')
end


%declare x_mode as global so that we start subsequent optimisations from
%the previous mode (speeds up nested optimisation).
global x_mode;
if isempty(x_mode)
  %no existing mode, use zeros as start
  x_mode = zeros(size(Qall,1),1);
end

%compute reorder (if needed)
p = amd(Qall + Aall'*Aall); %sparsity and reorder of Q_xy
Qall = Qall(p,p); %reorder
Aall = Aall(:,p);
x_mode = x_mode(p);

%find mode - using Newton-Raphson optimisation
%options = optimset('MaxIter', 1000);
x_mode = fminNR(@(x) GMRF_taylor_Po(x, y, Aall, Qall, E), x_mode);

%find the Laplace approximation of the denominator computed at the mode
[logp, ~, Q_xy] = GMRF_taylor_Po(x_mode, y, Aall, Qall, E);
%note that logp = -log_obs + x_mode'*Q*x_mode/2.

%Compute choleskey factor of Q_xy
[R_xy, p_xy] = chol(Q_xy);
if p_xy~=0
  %choleskey factor fail -> (almost) semidefinite matrix -> 
  %-> det(Q) ~ 0 -> log(det(Q)) ~ -inf -> negloglike ~ inf
  %Set negloglike to a REALLY big value
  negloglike = realmax;
  return;
end

%note that logp = -log_obs + x_mode'*Q*x_mode/2.
negloglike = logp + sum(log(diag(R_xy))) - 0.5*Qall_det;

%inverse reorder before returning
x_mode(p) = x_mode;

%print diagnostic information (progress)
fprintf(1, 'fval: %11.4e\n', negloglike);
