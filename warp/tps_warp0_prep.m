function prep0=tps_warp0_prep(varargin)
% TSP_WARP0_PREP  Precompute data for use in  tps_warp0
%
%  prep0=tps_warp0_prep(p0,sz0,axis0,lambda)
%  prep0=tps_warp0_prep(prep,p1,sz1,axis1)
%  prep0=tsp_warp0_prep(prep,x0,y0)
%
%  In the second and third form,  prep  is a struct given by a previous
%  invocation of  tps_warp_prep.  All three calls must be used before
%  using prep0 as input to tps_warp0.
%
%  p0 = The zero-energy control point configuration, np-by-2 matrix.
%  lambda = smoothing parameter.
%           lambda=0 gives an interpolating spline.
%  x0, y0 = coordinate matrices or vectors.

% $Id: tps_warp0_prep.m 3074 2006-10-26 11:17:46Z finn $

if (~isstruct(varargin{1}))
  p0 = varargin{1};
  sz0 = varargin{2};
  if (nargin<3), axis0 = []; else axis0 = varargin{3}; end
  if (nargin<4), lambda = []; else lambda = varargin{4}; end
  if isempty(axis0), axis0 = [1 sz0(2) 1 sz0(1)]; end
  if isempty(lambda), lambda = 0; end
  prep0.p0 = p0;
  prep0.sz0 = sz0;
  prep0.axis0 = axis0;
  prep0.lambda = lambda;

  np = size(p0,1);
  prep0.np = np;
  
  r2 = zeros(np,np);
  p0x = p0(:,ones(1,np));
  p0y = p0(:,1+ones(1,np));
  r2 = (p0x-p0x').^2 + (p0y-p0y').^2;
  prep0.K = r2.*log(r2+(r2==0))/2;
  prep0.P = [ones(np,1), p0];
  prep0.L = [prep0.K+eye(np)*lambda, prep0.P; prep0.P', zeros(3,3)];
  prep0.iL = inv(prep0.L);
  
  return;
end

prep0 = varargin{1};

% Try to determine which call was made:
if (nargin>3), calc_p1 = 1;
elseif any(size(varargin{2}) ~= size(prep0.p0)), calc_p1 = 0;
elseif (size(varargin{3},1) ~= 1), calc_p1 = 0;
elseif any(size(varargin{2}) ~= size(varargin{3})), calc_p1 = 1;
else, calc_p1 = 1;
end

if (calc_p1)
  p1 = varargin{2};
  sz1 = varargin{3};
  if (nargin<4), axis1 = []; else axis1 = varargin{4}; end
  if isempty(axis1), axis1 = [1 sz1(2) 1 sz1(1)]; end
  prep0.p1 = p1;
  prep0.sz1 = sz1;
  prep0.axis1 = axis1;
else
  x0 = varargin{2};
  y0 = varargin{3};

  prep0.x0 = x0;
  prep0.y0 = y0;
  prep0.szxy0 = size(x0);
  n = prep0.szxy0(1)*prep0.szxy0(2);

  x0 = x0(:);
  y0 = y0(:);

  r2 = (x0(:,ones(1,prep0.np))-ones(n,1)*prep0.p0(:,1)').^2 +...
       (y0(:,ones(1,prep0.np))-ones(n,1)*prep0.p0(:,2)').^2;

  prep0.U = r2.*log(r2+(r2==0))/2;
  prep0.B = [prep0.U,ones(n,1),x0,y0]/prep0.L;
  prep0.B(:,end-2:end) = [];
  prep0.Bim = reshape(prep0.B,[prep0.szxy0,prep0.np]);
end
