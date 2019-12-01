function [x,param,Fisher]=gmrf_snake(y,x0,param)
% GMRF_SNAKE Estimate a closed shape using a GMRF-snake model.
%
%  x=gmrf_snake(y,n,param)
%    Specifies the number of landmarks to use, n
%  x=gmrf_snake(y,x0,param)
%    Specifies an initial landmark configuration, x0
%  [x,param,Fisher]=gmrf_snake(...)
%    Also returns the complete parameter structure used in the function
%    as well as the Fisher information, i.e. the second derivative matrix
%    of the negative log-likelihood function.
%
%  GMRF-snake model:
%    Expectation: m
%    Precision: Q = alpha0*A + alpha1*W1'*A1*W1 + alpha2*W2'*A*W2
%
%  y:     The image, can be multi-dimensional
%  n/x0:  The number of landmarks (n) or initial landmark locations (x0)
%  param: A struct collecting all additional parameters.
%         Include only the parameters you wish to change from their
%         default values (see example below).
%    z01: The log-likelihood-ratio-image, log(p0(y)./p1(y)).
%         If z01 is empty or not present, the G-parameters are used
%         to calculate z01 instead.
%
%    G: Parameters for Gaussian pixel distributions
%      G.mu0, G.mu1: The pixel expectations for outside/inside
%      G.Sigma0, G.Sigma1: The pixel covariances for outside/inside
%
%    m : The expected snake shape landmarks (default=zeros)
%
%    alpha: The snake parameters, [alpha0, alpha1, alpha2]
%      alpha0>0 pulls the snake to m
%      alpha1>0 pulls the snake derivatives to the m-derivatives
%      alpha2>0 pulls the snake 2:nd derivatives to the m-2:nd-derivatives
%      Units: "one over (pixels squared)",
%             i.e. 1/("std.dev. along the whole curve")^2
%      Default=[0,0,1/300^2]
%
%    method: 'simple'=Simple local directions
%            'separation'=Separate the directions into the local normal and
%                         tangent directions
%            Default: 'simple'
%
%    scale: 0:  Only estimate affine shape deformations, i.e. scaling, and
%               translation, and skewing (includes rotations).
%           >1: Estimate deformations only at the specified scale.
%               scale=10 allows moving every 10:th landmark, dragging
%               the in-between-landmarks along in a smooth fashion.
%           Default: 0
%
%    smooth: Controls the amount of smoothing in the derivative
%            calculations.  Default=2
%
%    maxiter: The maximum number of interations in the optimisation.
%             Default=100
%
%    tol: The tolerance for determining convergence, in pixels.
%         Default=0.01
%
%    plotflag: If ==0, no plotting will be done
%              If >=1, plots the snake and z01-image (default)
%              If >=2, also plots the line-search data
%              If >=3, also plots samples from the posterior
%                      distribution, using the Fisher information
%
%  Example:
%    load lab7data
%    y = img{1};
%    imagesc(y), axis xy, colormap gray
%    param.G.mu0 = 0; param.G.Sigma0 = 0.01;
%    param.G.mu1 = 1; param.G.Sigma1 = 0.01;
%    param.scale = 0;  x=gmrf_snake(y,100,param);
%    param.scale = 10; x=gmrf_snake(y,x,param);
%    param.scale = 5;  x=gmrf_snake(y,x,param);
%    param.scale = 1;  x=gmrf_snake(y,x,param);

% $Id: gmrf_snake.m 4586 2012-10-08 16:18:33Z johanl $

if (nargin<3), param = []; end

sz = [size(y,1),size(y,2)];
D = size(y,3);

if (length(x0)==1) % x0=the number of landmarks
  n = x0;
  x0 = [cos(2*pi*((0:n-1)'/n))*sz(2)*0.3+sz(2)*0.5,...
        sin(2*pi*((0:n-1)'/n))*sz(1)*0.3+sz(1)*0.5];
else
  n = size(x0,1);
end

param0.z01 = [];
param0.G = [];
param0.alpha = [0,0,1/300^2];
param0.m = zeros(n,2);
param0.method = 'simple';
param0.scale = 0;
param0.smooth = 2;
param0.maxiter = 100;
param0.tol = 0.01;
param0.plotflag = 1;

obsolete_syntax = (~isfield(param,'G') & isfield(param,'mu0'));

param = defstruct(param,param0);

if (obsolete_syntax)
  warning(sprintf(...
      ['You are using an obsolete syntax for the parameters;\n',...
       'Please use the form  param.G.mu0, param.G.Sigma0, etc.,\n',...
       'instead of the old form  param.mu0, param.Sigma0, etc.\n',...
       'I will accept your parameters now, but I may stop doing\n',...
       'that in the future.  /gmrf_snake']));
  param.G.mu0 = param.mu0;
  param.G.Sigma0 = param.Sigma0;
  param.G.mu1 = param.mu1;
  param.G.Sigma1 = param.Sigma1;
elseif isempty(param.z01) & isempty(param.G)
  error(sprintf(...
      ['You must supply either the entire z01-log-likelihood-ratio\n',...
       'in param.z01 or parameters for a Gaussian model, in param.G']))
end

if (param.scale>=floor(n/4)) % Avoid stability problems.
  error(sprintf(...
      ['The scale parameter (%i) must be less than one fourth of ',...
       'the number of landmarks (%i).'],param.scale,n));
end

A = speye(n)/n;
W1 = sparse([1:n,1:n],...
            [1:n,2:n,1],...
            [-ones(1,n),ones(1,n)]*n,...
            n,n);
W2 = sparse([1:n,1:n,1:n],...
            [n,1:n-1,1:n,2:n,1],...
            [ones(1,n),-2*ones(1,n),ones(1,n)]*n^2,...
            n,n);
param0.Q_ = (param.alpha(1)*A + ...
             param.alpha(2)*W1'*A*W1 + ...
             param.alpha(3)*W2'*A*W2);
param0.Q = [param0.Q_, sparse(n,n); sparse(n,n), param0.Q_];

param = defstruct(param,param0);
param.M = vec(param.m);


Y = colstack(y);
if (isempty(param.z01))
  tmp = (Y-ones(length(Y),1)*param.G.mu0);
  z0 = icolstack(sum(tmp.*(tmp/param.G.Sigma0),2),sz)/2+...
       (log(det(param.G.Sigma0))+D*log(2*pi))/2;
  tmp = (Y-ones(length(Y),1)*param.G.mu1);
  z1 = icolstack(sum(tmp.*(tmp/param.G.Sigma1),2),sz)/2+...
       (log(det(param.G.Sigma1))+D*log(2*pi))/2;
  param.z01 = z1-z0;
end

dw0 = [1 2 1]/4;
dd0 = conv2([-1 1],[1 1]/2);
dw = dw0;
dd = dd0;
for k=2:param.smooth
  dw = conv2(dw,dw0);
  dd = conv2(dd,dw0);
end

g0 = conv2(dw,dw');
g1 = conv2(dd,dw');
g2 = conv2(dw,dd');

G0 = conv2(param.z01,rot90(g0,2),'same');
G1 = conv2(param.z01,rot90(g1,2),'same');
G2 = conv2(param.z01,rot90(g2,2),'same');
G0(1:param.smooth,:) = ones(param.smooth,1)*G0(param.smooth+1,:);
G0(end-param.smooth+1:end,:) = ones(param.smooth,1)*G0(end-param.smooth,:);
G0(:,1:param.smooth) = G0(:,param.smooth+1)*ones(1,param.smooth);
G0(:,end-param.smooth+1:end) = G0(:,end-param.smooth)*ones(1,param.smooth);
G2(1:param.smooth,:) = ones(param.smooth,1)*G2(param.smooth+1,:);
G2(end-param.smooth+1:end,:) = ones(param.smooth,1)*G2(end-param.smooth,:);
G1(:,1:param.smooth) = G1(:,param.smooth+1)*ones(1,param.smooth);
G1(:,end-param.smooth+1:end) = G1(:,end-param.smooth)*ones(1,param.smooth);

x = x0;

if (param.plotflag>=1)
  subplot(1,max(param.plotflag,1),1)
  imagesc(param.z01)
  axis xy
  hold on
  plot(x(:,1),x(:,2),'-r')
  hold off
  title(sprintf('Iteration %i',0))
  drawnow
end

conj_grad.alpha = 0;
conj_grad.s = 0;

for loop=1:param.maxiter
  x_old = x;
  
  no = (x([2:n,1],:)-x([n,1:n-1],:))*[0 -1;1 0];
  no = no./(sqrt(sum(no.^2,2))*[1,1]);
  L = sqrt(sum((x([2:n,1],:)-x).^2,2));
  L = (L+L([n,1:n-1]))/2;
  X = vec(x);
  N = vec(no);
  
  B0 = interp2(1:sz(2),1:sz(1),G0,x(:,1),x(:,2));
  B1 = interp2(1:sz(2),1:sz(1),G1,x(:,1),x(:,2));
  B2 = interp2(1:sz(2),1:sz(1),G2,x(:,1),x(:,2));
  dfdX = param.Q*(X-param.M) + N.*[L.*B0;L.*B0];
  
  tmp = L.*(no(:,1).*B1+no(:,2).*B2);
  d2fdX2 = param.Q + ...
           [spdiags(tmp.*no(:,1).*no(:,1),0,n,n),...
            spdiags(tmp.*no(:,1).*no(:,2),0,n,n);...
            spdiags(tmp.*no(:,2).*no(:,1),0,n,n),...
            spdiags(tmp.*no(:,2).*no(:,2),0,n,n)];

  if (param.scale==0) % Global positioning, scaling, and skew
    E = zeros(2*n,6);
    % Scaling, two directions:
    E(1:n,1) = x(:,1)-mean(x(:,1));
    E((1:n)+n,2) = x(:,2)-mean(x(:,2));
    % Positioning, two directions:
    E(1:n,3) = 1;
    E((1:n)+n,4) = 1;
    % Skew, two directions:
    E(1:n,5) = E((1:n)+n,1);
    E((1:n)+n,5) = E(1:n,1);
    E(1:n,6) = E((1:n)+n,2);
    E((1:n)+n,6) = E(1:n,2);
    if (0==1)
      % Better to include skew than diagonal scaling.
      % Including both results in rank deficiency.
      % Diagonal scaling, both diagonals:
      E(1:n,7) = E(1:n,1)+E((1:n)+n,2);
      E((1:n)+n,7) = E(1:n,1)+E((1:n)+n,2);
      E(1:n,8) = E(1:n,1)-E((1:n)+n,2);
      E((1:n)+n,8) = -E(1:n,1)+E((1:n)+n,2);
      % Remove components to avoid rank deficiency:
      E(:,[6,8]) = [];
    end
  elseif (param.scale==1)  % Local modifications
    E = speye(2*n);
  else % Medium-range modifications
    switch (param.method)
     case 'simple',
      num = max(1,floor(n/param.scale));
      width = param.scale*2-1;
      m0 = ceil(rand(1,1)*width);
      E = sparse(2*n,2*num);
      for k=1:num
        m = mod((k-1)*param.scale+m0-1,n)+1;
        fix = logical(ones(n,1));
        fix(mod([m-width:m-1,...
                 m+1:m+width]-1,n)+1) = 0;
        Ek = zeros(n,1);
        Ek(m) = 1;
        Ek(~fix) = -param.Q_(~fix,~fix)\(param.Q_(~fix,fix)*Ek(fix));
        E(1:n,k) = Ek;
        E((1:n)+n,k+num) = Ek;
      end
     case 'separation',
      num = max(1,floor(n/param.scale));
      width = param.scale*2-1;
      m0 = ceil(rand(1,1)*width);
      E = sparse(2*n,2*num);
      for k=1:num
        m = mod((k-1)*param.scale+m0-1,n)+1;
        fix = logical(ones(n,1));
        fix(mod([m-width:m-1,...
                 m+1:m+width]-1,n)+1) = 0;
        Ek = zeros(n,1);
        Ek(m) = 1;
        Ek(~fix) = -param.Q_(~fix,~fix)\(param.Q_(~fix,fix)*Ek(fix));
        E(:,k) = vec(Ek*no(m,:));
        E(:,k+num) = ...
            [Ek;Ek].*...
            vec(x(mod((1:n)+1-1,n)+1,:) -...
                x(mod((1:n)-1-1,n)+1,:));
      end
     otherwise,
      error(sprintf('Unknown method: "%s"',param.method))
    end
  end

  dX = -E*((E'*d2fdX2*E)\(E'*dfdX));
  
  if (max(abs(dX))<=5) & (0)
    X = X + dX;
    t = 1;
    line_search = 0;
  else
    line_search = 1;
    if (1)
      dX = dX*10/max(abs(dX));
    else
      % Conj. grad.
      conj_grad.v = -E*E'*dfdX + conj_grad.alpha*conj_grad.s;
      dX = conj_grad.v;
    end
    
    tt = linspace(-1,1,41).^3;
    for ti=1:length(tt)
      t = tt(ti);
      
      X_ = X + t*dX;
      x_ = ivec(X_,2);
      x_(:,1) = max(1,min(sz(2),x_(:,1)));
      x_(:,2) = max(1,min(sz(1),x_(:,2)));
      
      tmp2 = conv2(indicshape(x_,sz),g0,'same');
      g(ti) = sum(sum(param.z01.*tmp2));
      f(ti) = 0.5*(X_-param.M)'*param.Q*(X_-param.M) + ...
              g(ti);
    end
    [tmp,ti] = min(f);
    t = tt(ti);
    if (1)
      X = X + t*dX;
    else % Conj. grad.
      conj_grad.s = t*conj_grad.v;
      X = X + t*dX;

      x = ivec(X,2);
      x(:,1) = max(1,min(sz(2),x(:,1)));
      x(:,2) = max(1,min(sz(1),x(:,2)));
      
      no = (x([2:n,1],:)-x([n,1:n-1],:))*[0 -1;1 0];
      no = no./(sqrt(sum(no.^2,2))*[1,1]);
      L = sqrt(sum((x([2:n,1],:)-x).^2,2));
      L = (L+L([n,1:n-1]))/2;
      X = vec(x);
      N = vec(no);
      
      B0 = interp2(1:sz(2),1:sz(1),G0,x(:,1),x(:,2));
      dfdX_new = param.Q*(X-param.M) + N.*[L.*B0;L.*B0];
      grad_diff = dfdX_new-dfdX;
      conj_grad.alpha = -(grad_diff'*dfdX_new)/(grad_diff'*conj_grad.s);
    end
  end

  x = ivec(X,2);
  x(:,1) = max(1,min(sz(2),x(:,1)));
  x(:,2) = max(1,min(sz(1),x(:,2)));
  
  no = (x([2:n,1],:)-x([n,1:n-1],:))*[0 -1;1 0];
  no = no./(sqrt(sum(no.^2,2))*[1,1]);
  L = sqrt(sum((x([2:n,1],:)-x).^2,2));
  L = (L+L([n,1:n-1]))/2;
  X = vec(x);
  N = vec(no);
  
  if (param.plotflag>=1)
    subplot(1,max(param.plotflag,1),1)
    imagesc(param.z01);
    axis xy
    hold on
    plot(x(:,1),x(:,2),'-r')
    hold off
    title(sprintf('Iteration %i',loop))
  end
  if (param.plotflag>=2) & (line_search)
    subplot(1,max(param.plotflag,2),2)
    plot(tt,f-min(f),...
         tt,tt*(dfdX'*dX)+0.5*tt.^2*(dX'*d2fdX2*dX));
    title(sprintf('t = %2.4f',t))
  end

  if (param.plotflag>=3)
    subplot(1,max(param.plotflag,3),3)
    imagesc(param.z01)
    axis xy
    if (min(eig(full(d2fdX2)))>0)
      R = chol(d2fdX2);
      x_sim = ivec(X*ones(1,10)+R\randn(2*n,10),2);
      hold on
      for k=1:size(x_sim,3)
        plot(x_sim(:,1,k),x_sim(:,2,k),'m')
      end
      hold off
    end
    title(sprintf('Average landmark standard deviation = %2.2f',...
                  sqrt(full(mean(diag(inv(d2fdX2)))))));
  end

  change(loop) = max(abs(x(:)-x_old(:)));

  if (param.plotflag>=1)
    drawnow
  end

  if (change(loop)<param.tol)
    break;
  end
end

Fisher = d2fdX2;
