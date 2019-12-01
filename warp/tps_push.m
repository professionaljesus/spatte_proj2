function [X1,prep]=tps_push(varargin)
% TPS_PUSH Computes a push warp
%
%  This function treats the control points  p0  of the base image  X0
%  as the zero-energy configuration, warping from X0 to X1.
%
%  [X1,prep] = tps_push(X0,sz1,p0,p1,axis0,axis1,lambda,subsample)
%  X1 = tps_push(prep,X0,sz1,axis1,subsample)
%
%  sz1: The 2D-size of X1.
%       Default: sz1 = [size(X0,1),size(X0,2)]
%  axis0,axis1: [xmin xmax ymin ymax] for the image axes
%               Default: axis0 = [0.5 sz0(2)+0.5 0.5 sz0(1)+0.5]
%                        axis1 = [0.5 sz1(2)+0.5 0.5 sz1(1)+0.5]
%  lambda: smoothing parameter.
%          lambda=0 gives an interpolating spline.
%          Default: lambda=0.
%  subsample: The subsampling resolution
%             Default: subsample=3
%
% SEE ALSO: tps_prep, tps_pull, fillholes

% $Id: tps_push.m 4378 2011-08-16 12:28:05Z johanl $

if isstruct(varargin{1})
  prep = varargin{1};
  X0 = varargin{2};
  sz0 = prep.sz0;
  p0 = prep.p0;
  p1 = prep.p1;
  axis0 = prep.axis0;
  lambda = prep.lambda;
  if (nargin<3), sz1 = []; else sz1 = varargin{3}; end
  if (nargin<4), axis1 = []; else axis1 = varargin{4}; end
  if (nargin<5), subsample = []; else subsample = varargin{5}; end
  if isempty(sz1), sz1 = sz0; end
  if isempty(axis1), axis1 = [1 sz1(2) 1 sz1(1)]+0.5*[-1,1,-1,1]; end
  if isempty(subsample), subsample = 3; end
else
  X0  = varargin{1};
  sz0 = [size(X0,1) size(X0,2)];
  sz1 = varargin{2};
  p0 = varargin{3};
  p1 = varargin{4};
  if (nargin<5), axis0 = []; else axis0 = varargin{5}; end
  if (nargin<6), axis1 = []; else axis1 = varargin{6}; end
  if (nargin<7), lambda = []; else lambda = varargin{7}; end
  if (nargin<8), subsample = []; else subsample = varargin{8}; end
  if isempty(sz1), sz1 = sz0; end
  if isempty(axis0), axis0 = [1 sz0(2) 1 sz0(1)]+0.5*[-1,1,-1,1]; end
  if isempty(axis1), axis1 = [1 sz1(2) 1 sz1(1)]+0.5*[-1,1,-1,1]; end
  if isempty(lambda), lambda = 0; end
  if isempty(subsample), subsample = 3; end
  prep = tps_prep(tps_prep(tps_prep(p0,lambda),sz0,axis0),p1);
end

n = prod(prep.sz0);
d = size(X0,3);

delta0 = [(axis0(2)-axis0(1))/sz0(2),(axis0(4)-axis0(3))/sz0(1)];

axis0_x = axis0(1:2)+[1,-1]*0.5*delta0(1);
axis0_y = axis0(3:4)+[1,-1]*0.5*delta0(2);
[x0_mid,y0_mid] = meshgrid(linspace(axis0_x(1),axis0_x(2),sz0(2)),...
                           linspace(axis0_y(1),axis0_y(2),sz0(1)));

X1 = zeros([sz1,d]);
weights = zeros(sz1);
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

    x1_int = round((x1_interp(:)-axis1(1))/(axis1(2)-axis1(1))*...
                   sz1(2)+0.5);
    y1_int = round((y1_interp(:)-axis1(3))/(axis1(4)-axis1(3))*...
                   sz1(1)+0.5);
    I = (x1_int-1)*sz1(1)+y1_int;
    J = ones(length(I),1);
    ok = (x1_int>=1) & (x1_int<=sz1(2)) & ...
         (y1_int>=1) & (y1_int<=sz1(1));
    
    for k=1:d
      X0_ = interp2(x0_mid,y0_mid,double(X0(:,:,k)),...
                    x0,y0,...
                    '*cubic');
      if (k==1)
        weights0_ = (~isnan(X0_))*1;
        weights0_ = weights0_(:);
        
      end
      X0_(isnan(X0_)) = 0;
      X0_ = X0_(:);
      
      %% Transform to X1:
      if (k==1)
        weights1_ = accumarray([I(ok),J(ok)],weights0_(ok),[sz1(1)*sz1(2),1]);
        weights = weights + reshape(weights1_,sz1);
      end
      X1_ = accumarray([I(ok),J(ok)],X0_(ok),[sz1(1)*sz1(2),1]);
      X1(:,:,k) = X1(:,:,k) + reshape(X1_,sz1);
    end
  end
end

for k=1:d
  X1_ = X1(:,:,k)./(weights+(weights==0));
  X1_(weights==0) = NaN;
  X1(:,:,k) = X1_;
end
