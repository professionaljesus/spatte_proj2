function [z,Mz,Mf,Mzf]=mrf_sim(z0,N,alpha,beta,iter,gt)
% MRF_SIM Simulate a samples for a MRF
%
%  The field x is represented by a mxnxK-matrix
%  z_ik = x_i==k
%
%  The model is
%    P( x(i)=k | x(j), j \in N_i ) =
%      exp( alpha_ik + beta_k * f_ik ) / (sum_k exp( ... ))
%  which is equivalent to
%    E( z_ik | f_ik ) =
%      exp( alpha_ik + beta_k f_ik ) ...
%      / (sum_k' exp( alpha_ik' + beta_k f_ik' ))
%  where
%    f_ik = #{neighbours=k}
%
%  z=mrf_sim(z0,N,alpha,beta,iter)
%
%  z0    mxnxK, the initial indicator image
%  N     (2a+1)x(2b+1), neighbour-pattern. (0/1)
%  alpha 1x1 or 1xK or mxnxK, the expectation-forcing parameters
%  beta  1x1 or 1xK or KxK, the depency parameter(s)
%  iter  1x1, number of iterations
%  gt    mxn, ground truth data (1..K indicating apriori known class for a
%        pixel, 0 indicating unknown class). NOT IMPLEMENTED! 
%
%  The call
%    [x,Mz,Mf,Mzf]=mrf_sim(...)
%  also produces the estimated expectations
%    Mz_ik  = E(z_ik)      approx mean( E(z_ik | f_ik) )
%    Mf_ik  = E(f_ik)      approx mean( f_ik )
%    Mzf_ik = E(z_ik f_ik) approx mean( f_ik E(z_ik | f_ik) )
%  where the means are taken over the iterations, for iter > 0.
%  For iter==0, x=startx, and the expectations are calculated for this x.
%
% Example:
%   N = [0 1 0;1 0 1;0 1 0];
%   z = mrf_sim(zeros(100,120,3),N,log([0.2 0.3 0.5]),0.5,100);
%   for iter=1:100
%     z = mrf_sim(z,N,log([0.2 0.2 0.6]),[0.9 0.9 0],1);
%     image(rgbimage(z))
%     drawnow
%   end

% $Id: mrf_sim.m 5289 2019-01-22 10:49:35Z johanl $

if nargin<6, gt=[]; end

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

alpha = make_K_im(alpha,[m,n,K]);
if (length(beta)==1)
  beta = eye(K)*beta;
elseif ((size(beta,1)==1) || (size(beta,2)==1))
  beta = diag(beta);
end

if (iter==0)
  z = z0;
  [Mz,f] = calc_Mz_cond(z,N,alpha,beta,1:m,1:n,gt);
  Mf = f;
  Mzf = Mz.*f;
  return
end

if (nargout>1)
  Mz = zeros(m,n,K);
  Mf = zeros(m,n,K);
  Mzf = zeros(m,n,K);
end

KK = repmat(reshape(1:K,[1,1,K]),[m,n,1]);

z = z0;
ij_ = [kron(1:a,ones(1,b)); kron(ones(1,a),1:b)];
for loop=1:iter
  for ij=ij_(:,randperm(a*b))
    I = ij(1):a:m;
    J = ij(2):b:n;
    [Mz_cond,f] = calc_Mz_cond(z,N,alpha,beta,I,J,gt);
    if (nargout>1)
      Mz(I,J,:) = Mz(I,J,:)+Mz_cond;
      %      Mf(I,J,:) = Mf(I,J,:)+f(I,J,:);
      Mzf(I,J,:) = Mzf(I,J,:)+Mz_cond.*f(I,J,:);
    end
    e = rand(length(I),length(J));
    
    x = 1 + sum( cumsum( Mz_cond, 3 ) < ...
                 repmat( e, [1,1,K] ), ...
                 3 );
    z(I,J,:) = (repmat(x,[1,1,K])==KK(I,J,:)).*1;
  end
end
if (nargout>1)
  Mz = Mz/iter;
  %  Mf = Mf/iter;
  for k=1:K, Mf(:,:,k) = conv2(Mz(:,:,k),N,'same'); end
  Mzf = Mzf/iter;
end


function [Mz_cond,f] = calc_Mz_cond(z,N,alpha,beta,I,J,gt)
% Mz_cond  = E(z_ik | f_ik) )

K = size(z,3);
f = z;
for k=1:K, f(:,:,k) = conv2(z(:,:,k),N,'same'); end
Mz_cond = exp(alpha(I,J,:)+...
              reshape(reshape(f(I,J,:),...
                              [length(I)*length(J),K])*beta,...
                      [length(I),length(J),K]));
Mz_cond = Mz_cond./repmat(sum(Mz_cond,3),[1,1,size(z,3)]);

if ~isempty(gt)
  %NOT IMPLEMENTED (see Home Assignment)
  warning('fms150:gtNotImplementedInmrf_sim',...
    ['Ground truth available, but no implementation is ',...
    'available in mrf_sim!'])
  %pick out gt for the current block of conditionally independent pixels
  gt = gt(I,J);
  %columnstack the two images (easier to manipulate pixels)
  sz = [size(Mz_cond,1) size(Mz_cond,2)];
  Mz_cond = colstack(Mz_cond);
  gt = colstack(gt);
  %use gt to update class probabilities in Mz_cond
  %ADD CODE HERE!
  
  %inverse colstack back to suitable size.
  Mz_cond = icolstack(Mz_cond,sz);
end


function im=make_K_im(v,sz)
if (length(v)==1)
  im = v*ones(sz);
elseif ((size(v,1)==sz(1)) && (size(v,2)==sz(2)))
  if (size(v,3)==1)
    im = repmat(reshape(v,[sz(1:2),1]),[1,1,sz(3)]);
  else % v should already be the correct size; do nothing.
    im = v;
  end
else
  im = repmat(reshape(v,[1,1,sz(3)]),[sz(1:2),1]);
end
