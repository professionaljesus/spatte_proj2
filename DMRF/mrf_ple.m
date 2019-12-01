function [alpha,beta,pl]=mrf_ple(varargin)
% MRF_PLE Pseudo-likelihood estimation of MRF parameters
%
% See mrf_sim for the model description.
%
% Estimate (alpha, beta) for a given image z:
%   [alpha,beta]=mrf_ple(opt,N,z)
%   [alpha,beta]=mrf_ple(opt,N,z,I,J)
% Estimate (alpha, beta) within an EM algorithm:
%   [alpha,beta]=mrf_ple(opt,N,Mz,Mf,Mzf)
%   [alpha,beta]=mrf_ple(opt,N,Mz,Mf,Mzf,I,J)
%
% opt : [-1/0/1 -1/0/1 0/1]
%       opt(1)==1 : estimate alpha_k-vector.
%       opt(2)==1 : estimate beta_k-vector.
%       opt(1)==0 is equivalent to no alpha-parameters in the model.
%       opt(2)==0 estimates a common beta for all k.
%       opt(1)==-1 requires an extra argument to be placed last in the
%                  list, specifying a known alpha-value/vector; see example.
%       opt(2)==-1 requires an extra argument to be placed last in the
%                  list, specifying a known beta-value/vector; see
%                  example.
%       If both opt(1) and opt(2) are -1, the alpha-vector should be
%       placed before the beta-vector in the parameter list.
%       opt(3)==1 A mean zero Gaussian prior is used for beta.
%       opt(3)==0 No prior is used.
%
% I,J : index vectors for extracting the subgroup
%       of pixels that should be included in the pseudo-likelihood.
%
% Examples:
%   N = [0 1 0;1 0 1;0 1 0];
%   sz = [100,120,3];
%   z = mrf_sim(zeros(sz),N,log([0.2 0.3 0.5]),0.9,20);
%   alpha=mrf_ple([1,-1],N,z,1:2:sz(1),1:2:sz(2),0.9) % known beta=0.9
%   [alpha11,beta11]=mrf_ple([1,0],N,z,1:2:sz(1),1:2:sz(2))
%   [alpha12,beta12]=mrf_ple([1,0],N,z,1:2:sz(1),2:2:sz(2))
%   [alpha21,beta21]=mrf_ple([1,0],N,z,2:2:sz(1),1:2:sz(2))
%   [alpha22,beta22]=mrf_ple([1,0],N,z,2:2:sz(1),2:2:sz(2))
%   alpha_hat = (alpha11+alpha12+alpha21+alpha22)/4
%   a_hat = exp(alpha_hat)/sum(exp(alpha_hat))
%   beta_hat  = (beta11+beta12+beta21+beta22)/4

% $Id: mrf_ple.m 4837 2014-12-10 11:13:39Z johanl $

opt = varargin{1};
betaoffset = 1; 
if (opt(1)==-1),
  argoffset = 1;
  betaoffset = betaoffset+1;
else
  argoffset = 0;
end
if (opt(2)==-1), argoffset = argoffset+1; end
if (length(opt)<3), opt(3) = 1; end
useprior = (opt(3)>0);

N = varargin{2};
sz = size(varargin{3});
[m,n,K] = size(varargin{3});
if (nargin==3+argoffset) || (size(varargin{4},3)==1)
  Mz = varargin{3};
  Mf = Mz;
  for k=1:size(Mz,3), Mf(:,:,k) = conv2(Mz(:,:,k),N,'same'); end
  Mzf = Mz.*Mf;
  if (nargin>3+argoffset)
    I = varargin{4};
    J = varargin{5};
    argoffset = 5;
  else
    I = 1:m;
    J = 1:n;
    argoffset = 3;
  end
else
  Mz = varargin{3};
  Mf = varargin{4};
  Mzf = varargin{5};
  if (nargin>5+argoffset)
    I = varargin{6};
    J = varargin{7};
    argoffset = 7;
  else
    I = 1:m;
    J = 1:n;
    argoffset = 5;
  end
end
if (opt(1)<0), alpha = varargin{argoffset+1}; end
if (opt(2)<0), beta = varargin{argoffset+betaoffset}; end

Mz = Mz(I,J,:);
Mf = Mf(I,J,:);
Mzf = Mzf(I,J,:);
sz = size(Mz);
[m,n,K] = size(Mz);

iMz = reshape(sum(sum(Mz,1),2),[1,K]);
iMzf = sum(sum(Mzf,1),2);

if (opt(1)<0)
  if (length(alpha)==1)
    alpha = alpha*ones(1,K);
  else
    alpha = alpha(:)';
  end
end

if (opt(2)<0)
  if (length(beta)==1)
    beta = beta*ones(1,K);
  else
    beta = beta(:)';
  end
end

loss = @(a,b) ...
       -sum(log(a(:)).*iMz(:)) ...
       -sum(b(:).*iMzf(:)) ...
       +sum(sum(log(sum(make_K_im(a,sz).*...
                        exp(make_K_im(b,sz).*Mf),3)),1),2) ...
       +b(:)'*(-log(a(:)/sum(a(:))).*iMz(:).*b(:))/(10^2)*useprior;
       % A simple prior model, drawing beta towards 0.

fminopt = optimset('LargeScale','off','Display','off');
if (opt(1)>0)
  trsf = @(th) exp([th(1:K-1),-sum(th(1:K-1))]);
  if (opt(2)>0)
    th = [zeros(K-1,1);zeros(K,1)];
    th = fminsearch(@(th) loss(trsf(th'),...
                               th(K:2*K-1)'),th,fminopt);
    a = trsf(th'); a = a/sum(a);
    beta = th(K:2*K-1)';
  elseif (opt(2)==0)
    th = [zeros(K-1,1);0];
    th = fminsearch(@(th) loss(trsf(th'),...
                               th(K)*ones(1,K)),th,fminopt);
    a = trsf(th'); a = a/sum(a);
    beta = th(K)*ones(1,K);
  else
    if (1) % Simple iteration method.
      a = a_iter(beta,iMz,Mf);
    else
      th = [zeros(K-1,1)];
      th = fminsearch(@(th) loss(trsf(th'),...
                                 beta),th,fminopt);
      a = trsf(th'); a = a/sum(a);
    end
  end
  alpha = log(a);
elseif (opt(1)==0)
  if (opt(2)>0)
    a = ones(1,K)/K;
    th = fminsearch(@(th) loss(a,th'),zeros(K,1),fminopt);
    beta = th';
  elseif (opt(2)==0)
    a = ones(1,K)/K;
    th = fminsearch(@(th) loss(a,th*ones(1,K)),0,fminopt);
    beta = th*ones(1,K);
  else
    a = ones(1,K)/K;
  end
  alpha = log(a);
else
  a = exp(alpha);
  if (opt(2)>0)
    th = fminsearch(@(th) loss(a,th'),zeros(K,1),fminopt);
    beta = th';
  elseif (opt(2)==0)
    th = fminsearch(@(th) loss(a,th*ones(1,K)),0,fminopt);
    beta = th*ones(1,K);
  end
end

pl = -loss(a,beta);


function im=make_K_im(v,sz);
if (length(v)==1)
  im = v*ones(sz);
elseif ((size(v,1)==sz(1)) && (size(v,2)==sz(2)))
  im = repmat(reshape(v,[sz(1:2),1]),[1,1,sz(3)]);
else
  im = repmat(reshape(v,[1,1,sz(3)]),[sz(1:2),1]);
end

function a=a_iter(beta,iMz,Mf)
sz = size(Mf);
tmp0 = exp(make_K_im(beta,sz).*Mf);
a = (iMz(:)'+ones(1,sz(3))/sz(3))/(sz(1)*sz(2)+1);
for loop=1:10
  tmp1 = make_K_im(a,sz).*tmp0;
  tmp2 = sum(sum(tmp0./make_K_im(sum(tmp1,3),sz),1),2);
  a = iMz(:)'./tmp2(:)';
  a = a/sum(a);
end


% Some documentation:
%
% lnPL = sum_k alpha_k iMz(k) + sum_k beta_k iMzf(k)
%         - sum_i log sum_k exp( alpha_k + beta_k Mf(i,k) )
%      = sum_k log(a_k) iMz(k) + sum_k beta_k iMzf(k)
%         - sum_i log sum_k a_k exp( beta_k Mf(i,k) )
% dlnPL/da_k0 = iMz(k0)/a_k0
%         - sum_i exp( beta_k0 Mf(i,k0)
%                 / sum_k a_k exp( beta_k Mf(i,k) )
% dlnPL/dbeta_k0 = iMzf(k0)
%         - sum_i a_k Mf(i,k0) exp( beta_k0 Mf(i,k0)
%                 / sum_k a_k exp( beta_k Mf(i,k) )
