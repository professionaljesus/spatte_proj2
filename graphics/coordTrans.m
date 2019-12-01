function [loc2D,loc3D,lambert] = coordTrans(loc,intype,outtype)
% COORDTRANS Transforms coordinates between long/lat, 3D and lambert.
%
%  Call as
%   loc_out = coordTrans(loc,intype,outtype)
%  or
%   [loc2D,loc3D,lambert] = coordTrans(loc,intype,outtype)
%
%  Possible values for intype and outtype
%   '2d', 'longlat', 'long', 'lat', '3d', 'lambert'
%
%  2D coordinates are assumed to be lat/long in either degrees or radians.
%   They will be transforemd to [-pi,pi]X[-pi/2,pi/2] by scaling and
%   translation (if the x-coordinates is [0,2*pi]).
%
%  longlat, long, lat are all synonyms for 2D.
%
%  3D coordiantes are assumed to be on a zero-centred unit sphere.
%
%  lambert are similair to 2D but with a asin-transform of the
%   y-coordinate (ensures equal area projection).

% $Id: coordTrans.m 4837 2014-12-10 11:13:39Z johanl $

%determine intype
if strcmpi(intype,'2d') || strcmpi(intype,'longlat') || ...
      strcmpi(intype,'long') || strcmpi(intype,'lat')
  intype = 2;
elseif strcmpi(intype,'3d')
  intype = 3;
elseif strcmpi(intype,'lambert')
  intype = 0;
else
  error('Unknown intype');
end

%allocate space
loc2D = zeros(size(loc,1),2);
loc3D = zeros(size(loc,1),3);
lambert = zeros(size(loc,1),2);

if intype==0 || intype==2
  %first transform to radians
  if any(max(loc)>[2*pi pi/2]) || any(min(loc)<[-pi -pi/2])
    if intype==0
      loc(:,1) = loc(:,1)*pi/180;
    else
      loc = loc*pi/180;
    end
    if any(max(loc)>[2*pi pi/2]) || any(min(loc)<[-pi -pi/2])
      error('Too large coordinates')
    end
  end
  %centre to [-pi,pi]*[-1,1] or [-pi,pi]*[-pi/2,pi/2]
  if max(loc(:,1))>pi
    loc(:,1) = loc(:,1)-pi;
  end
end

if intype==2 %2D-coords in
  %2D -> %2D
  loc2D = loc;
  %2D -> %3D
  loc3D(:,3) = sin(loc2D(:,2));
elseif intype==3 %3D-coords in
  %First ensure that coordinates scale to one
  if max(abs( sum(loc.^2,2)-1 )) > 1e-12
    %standardise
    loc = loc/mean( sqrt(sum((10*loc).^2,2)) );
  end
  %3D -> %3D
  loc3D = loc;
  %3D -> %2D
  loc2D(:,1) = atan2(loc3D(:,2),loc3D(:,1));
  loc2D(:,2) = asin(loc3D(:,3));
elseif intype==0 %lambert-coords in
  %lambert -> lambert
  lambert = loc;
  %lambert -> %2D
  loc2D = zeros(size(lambert,1),2);
  loc2D(:,1) = lambert(:,1);
  loc2D(:,2) = asin(lambert(:,2));
  %lambert -> %3D
  loc3D = zeros(size(loc,1),3);
  loc3D(:,3) = lambert(:,2);
end

if intype==2 || intype==3
  %common 2D&3D -> lambert
  lambert(:,1) = loc2D(:,1);
  lambert(:,2) = loc3D(:,3);
end

if intype==0 || intype==2
  %common 2D&lambert -> 3D
  %loc3D(:,1) = cos(loc2D(:,2))
  %or using pythagoras:
  loc3D(:,1) = sqrt(1-loc3D(:,3).^2);
  loc3D(:,2) = sin(loc2D(:,1)) .* loc3D(:,1);
  loc3D(:,1) = cos(loc2D(:,1)) .* loc3D(:,1);
end

if nargout~=1
  return
end
if strcmpi(outtype,'2d') || strcmpi(outtype,'longlat') || ...
      strcmpi(outtype,'long') || strcmpi(outtype,'lat')
  return
elseif strcmpi(outtype,'3d')
  loc2D = loc3D;
elseif strcmpi(outtype,'lambert')
  loc2D = lambert;
else
  error('Unknown outtype');
end
