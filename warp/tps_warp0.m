function X0=tps_warp0(varargin)
% TPS_WARP0 Deform an image using a TPS warp.
%
%  This function treats the control points  p0  of the target image  X0
%  as the zero-energy configuration, warping from X1 to X0.
%
%  X0 = tps_warp0(X1,p0,p1)
%  X0 = tps_warp0(X1,p0,p1,sz0,sz1)
%  X0 = tps_warp0(X1,p0,p1,sz0,sz1,axis0,axis1)
%  X0 = tps_warp0(X1,p0,p1,sz0,sz1,axis0,axis1,lambda)
%
%  sz0, sz1: The 2D-size of X0 and X1.
%            Default: sz0 = sz1 = [size(X1,1),size(X1,2)]
%  axis0, axis1: [xmin xmax ymin ymax] for the image axes
%                Default: axis0 = [1 sz0(2) 1 sz0(1)]
%                         axis1 = [1 sz1(2) 1 sz1(1)]
%  lambda: smoothing parameter.
%          lambda=0 gives an interpolating spline.
%          Default: lambda=0.
%
% SEE ALSO: tps_warp1, tps_warp0_prep, fillholes
%
% Advanced usage:
%  X0 = tps_warp0(X1,prep0)
%  where prep0 must have been calculated by a sequence of calls to
%   tps_warp0_prep
%
% Example:
%  [x0,y0] = meshgrid(1:size(X0,2),1:size(X0,1));
%  X0 = tps_warp0(X1,tps_warp0_prep(...
%                    tps_warp0_prep(...
%                    tps_warp0_prep(p0,size(X0)),...
%                                   p1,size(X1)),...
%                                   x0,y0);

% $Id: tps_warp0.m 2957 2006-09-25 07:20:23Z johanl $

X1 = varargin{1};
if isstruct(varargin{2})
  prep0 = varargin{2};
else
  X1 = varargin{1};
  p0 = varargin{2};
  p1 = varargin{3};
  if (nargin<4), sz0 = []; else sz0 = varargin{4}; end
  if (nargin<5), sz1 = []; else sz1 = varargin{5}; end
  if (nargin<6), axis0 = []; else axis0 = varargin{6}; end
  if (nargin<7), axis1 = []; else axis1 = varargin{7}; end
  if (nargin<8), lambda = []; else lambda = varargin{8}; end
  if isempty(sz1), sz1 = [size(X1,1) size(X1,2)]; end
  if isempty(sz0), sz0 = sz1; end
  if isempty(axis0), axis0 = [1 sz0(2) 1 sz0(1)]; end
  if isempty(axis1), axis1 = [1 sz1(2) 1 sz1(1)]; end
  if isempty(lambda), lambda = 0; end
  [x0,y0] = meshgrid(linspace(axis0(1),axis0(2),sz0(2)),...
                     linspace(axis0(3),axis0(4),sz0(1)));
  prep0 = tps_warp0_prep(...
          tps_warp0_prep(...
          tps_warp0_prep(p0,sz0,axis0),...
                         p1,sz1,axis1),...
                         x0,y0);
end

n = prod(prep0.szxy0);
d = size(X1,3);

xy1 = reshape(prep0.B*prep0.p1,[prep0.szxy0,2]);
[x1,y1] = meshgrid(linspace(prep0.axis1(1),prep0.axis1(2),prep0.sz1(2)),...
		   linspace(prep0.axis1(3),prep0.axis1(4),prep0.sz1(1)));

X0 = zeros([prep0.szxy0,d]);
for k=1:d
  X0(:,:,k) = interp2(x1,y1,double(X1(:,:,k)),...
		      xy1(:,:,1),...
		      xy1(:,:,2),...
		      'linear');
end
