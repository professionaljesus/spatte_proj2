function X1=tps_warp1(X0,prep1)
% TPS_WARP1 Deform an image using an inverse TPS warp.
%
%  This function treats the control points  p0  of the original image
%  X0, as the zero-energy configuration, warping from X0 to X1 by
%  inverting the TPS coordinate function.
%
%  X1 = tps_warp1(X0,prep1)
%  where prep1 must have been calculated by  tps_warp1_prep
%  See  help tps_warp1_prep  for further information.
%
%  Example:
%   X1 = tps_warp1(X0,tps_warp1_prep(...
%                  tsp_warp1_prep(p0,p1,size(X0),size(X1))));
%
% SEE ALSO: tps_warp0, tps_warp1_prep, fillholes

% $Id: tps_warp1.m 2957 2006-09-25 07:20:23Z johanl $

if (nargin<4)
  axis1 = [];
end

n = prod(prep1.prep0.szxy0);
d = size(X0,3);

xy0 = reshape(prep1.xy0,[prep1.prep0.sz1,2]);
[x0,y0] = meshgrid(linspace(prep1.prep0.axis0(1),...
                            prep1.prep0.axis0(2),...
                            prep1.prep0.sz0(2)),...
		   linspace(prep1.prep0.axis0(3),...
                            prep1.prep0.axis0(4),...
                            prep1.prep0.sz0(1)));

X1 = zeros([prep1.prep0.sz1,d]);
for k=1:d
  X1(:,:,k) = interp2(x0,y0,X0(:,:,k),...
		      xy0(:,:,1),...
		      xy0(:,:,2),...
		      'linear');
end
