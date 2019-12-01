function [negl, Dl, D2l] = mrf_negLogPL(alpha, beta, z, f, beta_prior)
% MRF_NEGLOGPL Computes the negative log pseudo-likelihood of a MRF.
%
%  [l, Dl, D2l] = mrf_negLogPL(alpha, beta, z, f/N, beta_prior)
%
%  alpha [] or 1x(K-1), the expectation-forcing parameters. (For
%        identifiability the first element of alpha is taken as 0
%        i.e. alpha = [0 alpha];
%  beta  1x1 or 1xK, the depency parameter(s)
%  z     mxnxK, the field as an indicator image
%  f     Number of neighbours
%  N     (2a+1)x(2b+1), neighbour-pattern. (0/1)
%  beta_prior KxK PRECISION matrix for the prior on beta,
%             beta~N(0,beta_prior), if given as a sclar then
%             beta_prior*eye(K) is used.
%             beta_prior = 1/10 (i.e. variance of 10 is a reasonable choice)
%
%  The field x is represented by a mxnxK-matrix
%   z_ik = x_i==k
%
%  The model is
%    P( x(i)=k | x(j), j \in N_i ) =
%      exp( alpha_k + beta_k * f_ik ) / (sum_k exp( ... ))
%  where
%    f_ik = #{neighbours=k}
%  can be obtained from
%    [~,~,f] = mrf_sim(z, N, alpha, beta,0)
%
%  Returns logPL, along with first and second derivatives wrt alpha, beta.

% $Id: mrf_negLogPL.m 4837 2014-12-10 11:13:39Z johanl $

%% Default parameters
if nargin<5 || isempty(beta_prior), beta_prior=0; end

%expand alpha
if isempty(alpha), alpha=0; else alpha = [0 alpha(:)']; end

%is z and f of the same size
if ~isequal(size(z), size(f))
  %is not equal assume that f is a neighbourhood
  if ndims(z)~=3
    error('if z is colstacked then we need size(z)==size(f)');
  end
  %compute number of neighbours
  [~,~,f] = mrf_sim(z, f, alpha, beta, 0);
end

%check input sizes
if ndims(z)==3
  %columnstack z
  N = size(z,1)*size(z,2);
  K = size(z,3);
  z = reshape(z, [N K]);
  f = reshape(f, [N K]);
else
  K = size(z,2);
end

%sizes of alpha, beta and expand to K-by-1 matrices
if length(alpha)==1
  oneAlpha = true;
  alpha = repmat(alpha,[K 1]);
else
  oneAlpha = false;
  alpha = alpha(:);
end
if length(beta)==1
  oneBeta = true;
  beta = repmat(beta,[K 1]);
  %we're applying the prior to each beta and then summming, divide by
  %number of betas if only one beta.
  beta_prior = beta_prior/K;
else
  oneBeta = false;
  beta = beta(:);
end
%expand beta_prior to a K-by-K matrix
if all(size(beta_prior)==[1 1])
  beta_prior = beta_prior*eye(K);
end
if ~all(size(beta_prior)==[K K])
  error('beta_prior should be 1x1 or KxK')
end
%number of parameters
Npars = (K-1)*~oneAlpha + K*~oneBeta + oneBeta;

%compute z*f
zf = z.*f;
%compute sums
Sz = sum(z);
Szf = sum(zf);
%compute exp(alpha_k + beta_k*f_ik)
exp_alpha_beta = bsxfun(@times, f, beta');
exp_alpha_beta = bsxfun(@plus, exp_alpha_beta, alpha');
exp_alpha_beta = exp(exp_alpha_beta);
sum_exp = sum(exp_alpha_beta,2);
sum_log_exp = sum(log(sum_exp));

%compute neg logPL
negl = -(Sz*alpha + Szf*beta - sum_log_exp) + (beta'*beta_prior*beta)/2;

%compute derivative neg logPL
if nargout>1
  Dl = zeros(Npars,1);

  %derivatives wrt to alpha
  if ~oneAlpha
    Dl(1:(K-1)) = -Sz(2:end) ...
      + sum(bsxfun(@rdivide, exp_alpha_beta(:,2:end), sum_exp));
    offset = K;
  else
    offset = 1;
  end
  %derivatives wrt to beta
  dl_beta = (beta_prior*beta)' - Szf ...
      + sum(bsxfun(@rdivide, f.*exp_alpha_beta, sum_exp));
  if ~oneBeta
    Dl(offset:end) = dl_beta;
  else
    Dl(end) = sum(dl_beta);
  end
end

%compute second derivative neg logPL
if nargout>2
  D2l = zeros(Npars,Npars);
  %compute sum(exp(alpha+beta))^2
  sum_exp_2 = sum_exp.^2;
  %also compute sum(exp(alpha(-k)+beta(-k)))^2
  sum_exp_but_k = zeros(size(exp_alpha_beta));
  for k=1:K
    Ind = true(1,K);
    Ind(k)= false;
    sum_exp_but_k(:,k) = sum(exp_alpha_beta(:,Ind),2);
  end
  
  %derivatives wrt to alpha
  if ~oneAlpha
    for k=2:K
      D2l(k-1,k-1) = sum(exp_alpha_beta(:,k) .* ...
        sum_exp_but_k(:,k) ./ sum_exp_2);
      for l=(k+1):K
        D2l(k-1,l-1) = -sum(exp_alpha_beta(:,k) .* ...
          exp_alpha_beta(:,l) ./ sum_exp_2);
        D2l(l-1,k-1) = D2l(k-1,l-1);
      end
    end
    %cross derivatives
    D2l_beta = zeros(K,K-1);
    for k=2:K
      for l=1:K
        tmp = exp_alpha_beta(:,k).*f(:,l);
        if k==l
          tmp = tmp.*sum_exp_but_k(:,l);
        else
          tmp = -tmp.*exp_alpha_beta(:,l);
        end
        D2l_beta(l,k-1) = sum(tmp./sum_exp_2);
      end
    end
    if oneBeta, D2l_beta=sum(D2l_beta); end
    D2l(K:end,1:(K-1)) = D2l_beta;
    D2l(1:(K-1), K:end) = D2l_beta';
    offset = K;
  else
    offset = 1;
  end
  
  %derivatives wrt to beta
  D2l_beta = zeros(K,K);
  for k=1:K
    D2l_beta(k,k) = sum(f(:,k).*exp_alpha_beta(:,k) .* ...
        sum_exp_but_k(:,k) ./ sum_exp_2);
    if ~oneBeta
      for l=(k+1):K
        D2l_beta(k,l) = -sum(f(:,k).*exp_alpha_beta(:,k) .* ...
        f(:,l).*exp_alpha_beta(:,l) ./ sum_exp_2);
        D2l_beta(l,k) = D2l_beta(k,l);
      end
    end
  end
  D2l_beta = D2l_beta + beta_prior;
  if oneBeta
    D2l(offset:end,offset:end) = sum(diag(D2l_beta));
  else
    D2l(offset:end,offset:end) = D2l_beta;
  end
end
