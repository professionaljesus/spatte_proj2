function globe_plot(x,points,map,globe)
% GLOBE_PLOT Plot a field defined on a spherical grid 
%
% globe_plot(x)
% globe_plot(x,[],[],[])
% globe_plot(x,points,map,globe)
%
%   x is an m-by-n matrix of field values, corresponding to locations
%     in a grid with equal degree distances between all rows and between
%     all columns.  The first row is the north pole, and row m is the
%     south pole.  The first column is at longitude -180 degrees.
%   If present, points.long and points.lat  contain the longitudes
%     and latitudes of points that should be marked on the plot.
%   If present, map.long and map.lat  contain the longitudes and
%     latitudes of coastlines that should be marked on the plot.
%   globe==0 gives a planar projection with equal degree differences
%            between all rows and between all columns.
%   globe==1 gives a spherical plot. (default)
%
% Examples:
% % Load project data and global map:
% load proj2.mat
% % Plot random fields:
% globe_plot(randn(30,60),[],globe_map,0)
% globe_plot(randn(30,60),[],globe_map,1)
% % Estimate the temperature field X:
% ...
% % Plot the result:
% x = icolstack(X,sz);
% globe_plot(x,loc,globe_map,0)
% globe_plot(x,loc,globe_map,1)
% axis equal
% axis vis3d
% rotate3d on

% $Id: globe_plot.m 4586 2012-10-08 16:18:33Z johanl $

if (nargin<4), globe = []; end
if (nargin<3), map = []; end
if (nargin<2), points = []; end

if isempty(globe), globe = 1; end
if isempty(map), map = []; end
if isempty(points), points = []; end

sz = size(x);

if (globe == 1)
  lon = linspace(-180,180,sz(2)+1)/360*2*pi;
  lat = linspace(90,-90,sz(1))/360*2*pi;

  [la,lo] = ndgrid(lat,lon);
  [Px,Py,Pz] = sph2cart(lo,la,1);
  
  surf(Px,Py,Pz,x(:,[1:end,1]))
  shading interp
  axis equal
  
  if (~isempty(map))
    [Px,Py,Pz] = sph2cart(map.long/360*2*pi,...
                          map.lat/360*2*pi,1);
    hold on
    plot3(Px*1.001,Py*1.001,Pz*1.001,'k')
    hold off
  end
  
  if (~isempty(points))
    [Px,Py,Pz] = sph2cart(points.long/360*2*pi,...
                          points.lat/360*2*pi,1);
    hold on
    plot3(Px*1.001,Py*1.001,Pz*1.001,'.m')
    hold off
  end
elseif (globe==0)
  lon = linspace(-180,180,sz(2)+1);
  lat = linspace(90,-90,sz(1));
  [la,lo] = ndgrid(lat,lon);

  surf(lo,la,lo*0,x(:,[1:end,1]))
  shading interp
  axis([-180,180,-90,90])
  
  if (~isempty(map))
    hold on
    plot3(map.long,map.lat,map.long*0+0.001,'k')
    hold off
  end
  
  if (~isempty(points))
    hold on
    plot3(points.long,points.lat,points.long*0+0.002,'.m')
    hold off
  end
else % globe==2, area-correct projection
  lon = linspace(-180,180,sz(2)+1);
  lat = linspace(90,-90,sz(1));
  [la,lo] = ndgrid(lat,lon);

  surf(lo,sin(la/360*2*pi),lo*0,x(:,[1:end,1]))
  shading interp
  axis([-180,180,-1,1])
  
  if (~isempty(map))
    hold on
    plot3(map.long,sin(map.lat/360*2*pi),map.long*0+0.001,'k')
    hold off
  end
  
  if (~isempty(points))
    hold on
    plot3(points.long,sin(points.lat/360*2*pi),points.long*0+0.001,'.m')
    hold off
  end
end
