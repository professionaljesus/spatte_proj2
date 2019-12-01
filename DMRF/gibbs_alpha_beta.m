function [alpha, beta, acc] = gibbs_alpha_beta(alpha, beta, z, f, beta_prior, MHsigma)
% GIBBS_ALPHA_BETA Samples a new set of alpha,beta given z using a MH-gibbs step
%
%  [alpha, beta, acc] = gibbs_alpha_beta(alpha, beta, z, f, beta_prior, MHsigma)
%  [alpha, beta, acc] = gibbs_alpha_beta(alpha, beta, z, N, beta_prior, MHsigma)
%
%  alpha:[] or 1x(K-1), the expectation-forcing parameters. For
%        identifiability the first element of alpha is taken as 0
%        i.e. the alpha used in mrf_sim should be extended as [0 alpha];
%  beta: 1x1 or 1xK, the depency parameter(s)
%  z:    mxnxK, the field as an indicator image
%  f:    Number of neighbours
%  N:    (2a+1)x(2b+1), neighbour-pattern. (0/1)
%  beta_prior: KxK PRECISION matrix for the prior on beta,
%              beta~N(0,beta_prior), if given as a sclar then
%              beta_prior*eye(K) is used.
%              beta_prior = 1/10 (i.e. variance of 10 is a reasonable choice)
%  MHsigma: proposal variance in the Metropolis Hastings step (defaults to 1e-4)
%
%  The field x is represented by a mxnxK-matrix
%   z_ik = x_i==k
%
%  The model is
%    P( x(i)=k | x(j), j \in N_i ) =
%      exp( alpha_k + beta_k * f_ik ) / (sum_k exp( ... ))
%  where
%    f_ik = #{neighbours=k}
%  can be obtained from (recall the need to extend alpha)
%    [~,~,f] = mrf_sim(z, N, [0 alpha], beta,0)
%
%  Returns updated alpha, beta and acceptance rate.

% $Id: gibbs_alpha_beta.m 5254 2018-12-09 18:23:49Z johanl $

%% Default parameters
if nargin<5 || isempty(beta_prior), beta_prior=[]; end
if nargin<6 || isempty(MHsigma), MHsigma=1e-4; end

%check input sizes
if ndims(z)==3
  K = size(z,3);
else
  K = size(z,2);
end

%sizes of alpha, beta and expand to K-by-1 matrices
oneAlpha = isempty(alpha);
oneBeta = (length(beta)==1);
%number of parameters
Npars = (K-1)*~oneAlpha + K*~oneBeta + oneBeta;

if numel(MHsigma)==1, MHsigma=eye(Npars)*MHsigma; end

%% Create neglogLike function to sample from
if oneAlpha
  logL = @(x) mrf_negLogPL([], x, z, f, beta_prior);
else
  logL = @(x) mrf_negLogPL(x(1:(K-1)), x(K:end), z, f, beta_prior);
end

%% Sample
[alpha_beta, acc] = HM_sampling(logL, [alpha(:);beta(:)], MHsigma);

if oneAlpha
  alpha = [];
  beta = alpha_beta;
else
  alpha = reshape(alpha_beta(1:(K-1)),1,[]);
  beta = reshape(alpha_beta(K:end),1,[]);
end
    
    
function [x, acc] = HM_sampling(neglogL, x, theta)
logL_old = -neglogL(x);

%sample new values
R = chol(theta*( eye(length(x)))); 
x_new = x + R'*randn(size(x));

%compute acceptance probability
logL_new = -neglogL(x_new);
alpha = min(exp(logL_new - logL_old), 1);
if rand(1) < alpha
   x = x_new;
   acc = 1;
else
  acc  = 0;
end
