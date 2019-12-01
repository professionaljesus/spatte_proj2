function [mx,vx,cx]=calc_gmrf_props(Q,b,c,method)
% CALC_GMRF_PROPS  Helper function for calculating pointwise GMRF properties.
%
%  mx = calc_gmrf_props(Q,b,c,method)
%  [mx,vx] = calc_gmrf_props(Q,b,c,method)
%  [mx,vx,cx] = calc_gmrf_props(Q,b,c,method)
%
%  Q = A precision matrix
%  b = a vector (See the canonical representation)
%  c = a sparse matrix
%      The actual matrices to supply depend on the circumstances!
%      In HA2 c = A'*A/sigma2_epsilon.
%  method = 1 for exact calculations using back-substitution (default)
%           2 for approximate computations using simulation
%           Which method is faster depends on the size of the problem.
%
%  mx = (Q+c)^{-1} b (used for expectations)
%  vx(i) = element (i,i) of (Q+c)^{-1}  (used for pointwise variances)
%  cx(i,j) = element (i,j) of (Q+c)^{-1} for (i,j) with Q(i,j)~=0
%            (not used in HA2)
%
%  For method=2 the variance and covariance (vx and cx) are approximated from
%  simulations and results will vary between calls.
%
%  Input and default values:
%    Q sparse precision matrix, nxn
%    b default = zeros(n,1)
%    c default = sparse(n,n)
%    method default = 1

n = size(Q,1);

if (nargin<2), b = []; end
if (nargin<3), c = []; end
if (nargin<4), method = []; end
if isempty(b), b = zeros(n,1); end
if isempty(c), c = sparse(n,n); end
if isempty(method), method = 1; end % default = Back-substitution

T0 = cputime;
R = chol(Q+c);
T = cputime;
mx = R\((b'/R)');

if (nargout<=1), return, end

switch (method)
  case 1, % Back-substitution
    vx = zeros(n,1);
    if (nargout>=3)
      cx = sparse([],[],[],n,n,nnz(Q+c));
    end
    for i=1:n
      e = zeros(1,n);
      e(i) = 1;
      eR = e/R;
      vx(i) = sum(eR.^2);
      if (nargout>=3)
        cx(i,i) = vx(i);
        j = find(Q(i,i+1:end)+c(i,i+1:end));
        f = sparse((1:length(j))',j,ones(length(j),1),length(j),n);
        cxs = ((f/R)*eR')';
        cx(i,j) = cxs;
        cx(j,i) = cxs;
      end
    end
  case 2, % Simulate
    k = 0;
    sum_x = zeros(n,1);
    sum_x2 = zeros(n,1);
    if (nargout>=3)
      [i,j] = find(Q);
      sum_xx = zeros(length(i),1);
    end
    tic
    while ((cputime-T)<(T-T0+0.01)*100)
      k = k+1;
      x = R\randn(n,1);
      sum_x2 = sum_x2+x.^2;
      if (nargout>=3)
        sum_xx = sum_xx+x(i).*x(j);
      end
    end
    
    vx = sum_x2/k;
    if (nargout>=3)
      cx = sparse(i,j,sum_xx/k,n,n);
    end
end
