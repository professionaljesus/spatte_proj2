function zhat=mrf_icm(z0,N,alpha,beta,maxiter)
% MRF_ICM Estimate the MAP field from an MRF model
%
%  See mrf_sim for the model description.
%
%  zhat=mrf_icm(z0,N,alpha)
%  zhat=mrf_icm(z0,N,alpha,beta)
%  zhat=mrf_icm(z0,N,alpha,beta,maxiter)
%
%  z0      mxnxK, an initial indicator image
%  N       (2a+1)x(2b+1), neighbour-pattern. (0/1)
%  alpha   mxnxK, the expectation-forcing parameters
%  beta    1x1 or 1xK or KxK, the depency parameter(s) (default=0)
%  maxiter 1x1, the maximal number of iterations (default=Inf)
%
% Examples:
%   alpha = log("data-likelihood");  % Calculate posterior alpha-values.
%   zhat = mrf_icm(z0,N,alpha,beta); % Calculate point estimate of z.
%   % View the optimisation:
%   zhat = z0;
%   for loop=1:100 
%     zhat = mrf_icm(zhat,N,alpha,beta,1); % maxiter==1
%     imagesc(rgbimage(zhat))
%     drawnow
%   end

% $Id: mrf_icm.m 4837 2014-12-10 11:13:39Z johanl $

if (nargin<4), beta = []; end
if isempty(beta), beta = 0; end
if (nargin<5), maxiter = []; end
if isempty(maxiter), maxiter = Inf; end

[m,n,K] = size(z0);
if sum(N-rot90(N,2))
  error('The neighbourhood must have reflective symmetry.')
end
[a,b] = size(N);
if ((mod(a,2)~=1) || (mod(b,2)~=1))
  error('The neighbourhood must have odd width and height.')
end
a = ceil(a/2);
b = ceil(b/2);

KK = repmat(reshape(1:K,[1,1,K]),[m,n,1]);

zprevious = z0*0;
zhat = z0;
loop = 0;
ij_ = [kron(1:a,ones(1,b)), kron(a:-1:1,ones(1,b));...
       kron(ones(1,a),1:b), kron(ones(1,a),b:-1:1)];
ij_(:,[a*b,2*a*b]) = [];
while (loop<maxiter) && any(zhat(:)~=zprevious(:))
  loop = loop+1;
  zprevious = zhat;
  for ij=ij_
    I = ij(1):a:m;
    J = ij(2):a:n;
    [unused,Mz] = mrf_sim(zhat,N,alpha,beta,0);
    [unused,xhat] = max(Mz(I,J,:),[],3);
    zhat(I,J,:) = (repmat(xhat,[1,1,K])==KK(I,J,:)).*1;    
  end  
end
