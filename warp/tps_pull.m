function [X0,prep]=tps_pull(varargin)
% TPS_PULL Computes a pull warp
%
%  This function treats the control points  p0  of the target image  X0
%  as the zero-energy configuration, warping from X1 to X0.
%
%  [X0,prep] = tps_pull(sz0,X1,p0,p1,axis0,axis1,lambda,subsample)
%  X0 = tps_pull(prep,X1,axis1,subsample)
%
%  sz0: The 2D-size of X0.
%       Default: sz0 = [size(X1,1),size(X1,2)]
%  axis0,axis1: [xmin xmax ymin ymax] for the image axes
%               Default: axis0 = [0.5 sz0(2)+0.5 0.5 sz0(1)+0.5]
%                        axis1 = [0.5 sz1(2)+0.5 0.5 sz1(1)+0.5]
%  lambda: smoothing parameter.
%          lambda=0 gives an interpolating spline.
%          Default: lambda=0.
%  subsample: The subsampling resolution
%             Default: subsample=3
%
% SEE ALSO: tps_prep, tps_push, fillholes

% $Id: tps_pull.m 4378 2011-08-16 12:28:05Z johanl $

if isstruct(varargin{1})
  prep = varargin{1};
  sz0 = prep.sz0;
  X1 = varargin{2};
  sz1 = [size(X1,1) size(X1,2)];
  p0 = prep.p0;
  p1 = prep.p1;
  axis0 = prep.axis0;
  lambda = prep.lambda;
  if (nargin<3), axis1 = []; else axis1 = varargin{3}; end
  if (nargin<4), subsample = []; else subsample = varargin{4}; end
  if isempty(axis1), axis1 = [1 sz1(2) 1 sz1(1)]+0.5*[-1,1,-1,1]; end
  if isempty(subsample), subsample = 3; end
else
  sz0 = varargin{1};
  X1  = varargin{2};
  p0 = varargin{3};
  p1 = varargin{4};
  if (nargin<5), axis0 = []; else axis0 = varargin{5}; end
  if (nargin<6), axis1 = []; else axis1 = varargin{6}; end
  if (nargin<7), lambda = []; else lambda = varargin{7}; end
  if (nargin<8), subsample = []; else subsample = varargin{8}; end
  sz1 = [size(X1,1) size(X1,2)];
  if isempty(sz0), sz0 = sz1; end
  if isempty(axis0), axis0 = [1 sz0(2) 1 sz0(1)]+0.5*[-1,1,-1,1]; end
  if isempty(axis1), axis1 = [1 sz1(2) 1 sz1(1)]+0.5*[-1,1,-1,1]; end
  if isempty(lambda), lambda = 0; end
  if isempty(subsample), subsample = 3; end
  prep = tps_prep(tps_prep(tps_prep(p0,lambda),sz0,axis0),p1);
end

n = prod(prep.sz0);
d = size(X1,3);

delta0 = [(axis0(2)-axis0(1))/sz0(2),(axis0(4)-axis0(3))/sz0(1)];

axis1_x = axis1(1:2)+[1,-1]*0.5*(axis1(2)-axis1(1))/sz1(2);
axis1_y = axis1(3:4)+[1,-1]*0.5*(axis1(4)-axis1(3))/sz1(1);
[x1,y1] = meshgrid(linspace(axis1_x(1),axis1_x(2),sz1(2)),...
		   linspace(axis1_y(1),axis1_y(2),sz1(1)));

X0 = zeros([prep.sz0,d]);
weights = zeros(prep.sz0);
xy0 = zeros([prep.sz0,2]);
axis0_x_base = prep.axis0(1:2)+[0,-1]*delta0(1);
axis0_y_base = prep.axis0(3:4)+[0,-1]*delta0(2);
for subsample_x=1:subsample
  axis0_x = axis0_x_base + (2*subsample_x-1)/(2*subsample)*delta0(1);
  for subsample_y=1:subsample
    axis0_y = axis0_y_base + (2*subsample_y-1)/(2*subsample)*delta0(2);
    [x0,y0] = meshgrid(linspace(axis0_x(1),axis0_x(2),prep.sz0(2)),...
                       linspace(axis0_y(1),axis0_y(2),prep.sz0(1)));
    xy0(:,:,1) = x0;
    xy0(:,:,2) = y0;
    prep_interp = tps_prep(prep,xy0);
    x1_interp = prep_interp.xy1(:,:,1);
    y1_interp = prep_interp.xy1(:,:,2);
    for k=1:d
      X0_ = interp2(x1,y1,double(X1(:,:,k)),...
                    x1_interp,y1_interp,...
                    '*cubic');
      if (k==1)
        weights = weights + (~isnan(X0_));
      end
      X0_(isnan(X0_)) = 0;
      X0(:,:,k) = X0(:,:,k) + X0_;
    end
  end
end

for k=1:d
  X0_ = X0(:,:,k)./(weights+(weights==0));
  X0_(weights==0) = NaN;
  X0(:,:,k) = X0_;
end
