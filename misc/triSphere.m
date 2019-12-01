function [Xn,T,a,A] = triSphere(X, angle_thres, area_thres, plotF)
% TRISPHERE Triangulates a sphere.
%
%  [Xn,T,a,A] = triSphere(X, angle_thres, area_thres, plotF)
%
%  X is a n-by-3 matrix of points on the sphere assumed to have radie 1.
%  angle_thres is the minimum angle of any triangle
%  area_thres is the maximum area of any triangle
%  plotF indicates if plots of the progression are wanted.
%
%  The output consists of
%  Xn a new list of points, the first n points are equal to the input
%  points, other points can be added to comply with the angle and/or area
%  requirments.
%  T is a list of triangles
%  a it the minimum angle for each triangle.
%  A is the area of each triangle.

% $Id: triSphere.m 4586 2012-10-08 16:18:33Z johanl $

if nargin<2, angle_thres=[]; end;
if nargin<3, area_thres=[]; end;
if nargin<4, plotF=[]; end;
  
if isempty(angle_thres), angle_thres=0; end;
if isempty(area_thres), area_thres=0; end;
if isempty(plotF), plotF=0; end;
if angle_thres>pi/4, angle_thres = angle_thres*pi/180; end;
if angle_thres>(33*pi/180)
  error('Threshold on angle >33. Alogrithm will probably not converge');
end

%initate vector of added points
Xn = zeros(0,3);
%indicator if we have added points
Addpoint = true;

while Addpoint
  %create temporary X vector
  Xt = [X;Xn];
  %calculate triangulation, angles, areas and circumcentres
  [T,a,A,C] = calcTristats(Xt);
  if plotF>0
    %do plots
    doPlot(T,Xt,a,A,angle_thres,area_thres)
  end

  Addpoint = false;
  %add new points starting with the smallest angle
  while ~isempty(T) && angle_thres>0 && min(a)<angle_thres
    [tmp,Ind] = min(a);
    %first add circumcentre
    Xn = [Xn;C(Ind,:)];
    Addpoint = true;
    %find all triangles adjacent to this
    I = any(T==T(Ind,1),2) | any(T==T(Ind,2),2) | any(T==T(Ind,3),2);
    %drop these triangles
    T(I,:) = []; a(I) = []; A(I) = []; C(I,:) = [];
  end
  while ~isempty(T) && area_thres>0 && max(A)>area_thres
    [tmp,Ind] = max(A);
    %first add circumcentre
    Xn = [Xn;C(Ind,:)];
    Addpoint = true;
    %find all triangles adjacent to this
    I = any(T==T(Ind,1),2) | any(T==T(Ind,2),2) | any(T==T(Ind,3),2);
    %drop these triangles
    T(I,:) = []; a(I) = []; A(I) = []; C(I,:) = [];
  end
end
Xn = [X;Xn];


function [T,a,A,C] = calcTristats(X)
%first find a triangulation, i.e. the convexhull.
T = convhulln(X);

%exctract vectors between triangle nodes
v = zeros(size(T,1),size(T,2),3);
v(:,:,1) = X(T(:,2),:)-X(T(:,1),:);
v(:,:,2) = X(T(:,3),:)-X(T(:,2),:);
v(:,:,3) = X(T(:,1),:)-X(T(:,3),:);

%calculate minimum angle for each triangle
a = zeros(size(T,1),3);
a(:,1) = -sum(v(:,:,1).*v(:,:,3),2);
a(:,2) = -sum(v(:,:,1).*v(:,:,2),2);
a(:,3) = -sum(v(:,:,2).*v(:,:,3),2);
a(:,1) = a(:,1)./sqrt(sum(v(:,:,1).^2,2).*sum(v(:,:,3).^2,2));
a(:,2) = a(:,2)./sqrt(sum(v(:,:,1).^2,2).*sum(v(:,:,2).^2,2));
a(:,3) = a(:,3)./sqrt(sum(v(:,:,2).^2,2).*sum(v(:,:,3).^2,2));
a = min(acos(a),[],2);

%calculate area of all the triangles
A = [-v(:,2,1).*v(:,3,3)+v(:,3,1).*v(:,2,3) ...
     -v(:,3,1).*v(:,1,3)+v(:,1,1).*v(:,3,3) ...
     -v(:,1,1).*v(:,2,3)+v(:,2,1).*v(:,1,3)];
A = 1/2*sqrt(sum(A.^2,2));

%calculate the circumcentre of each triangle
alpha = sum(v(:,:,2).^2,2).*sum(-v(:,:,1).*v(:,:,3),2)./(8*A.^2);
beta = sum(v(:,:,3).^2,2).*sum(-v(:,:,1).*v(:,:,2),2)./(8*A.^2);
gamma = sum(v(:,:,1).^2,2).*sum(-v(:,:,3).*v(:,:,2),2)./(8*A.^2);
alpha = repmat(alpha,[1 3]);
beta = repmat(beta,[1 3]);
gamma = repmat(gamma,[1 3]);
C = alpha.*X(T(:,1),:) + beta.*X(T(:,2),:) + gamma.*X(T(:,3),:);
C = C./repmat(sqrt(sum(C.^2,2)),[1 3]);


function doPlot(T,X,a,A,thres_a,thres_A)
%construct some figures
subplot(221)
trimesh(T,X(:,1),X(:,2),X(:,3),-inf(size(X,1),1));
axis equal
subplot(223)
hist(a*180/pi,50);
axis tight
AX = axis;
if thres_a>0
  hold on
  plot(thres_a*180/pi*[1 1],AX([3 4]),'r')
  hold off
end
axis([0 AX(2:4)])
title( min(a)*180/pi );

subplot(224)
hist(log10(A),50);
axis tight
AX = axis;
if thres_A>0
  hold on
  plot(log10(thres_A)*[1 1],AX([3 4]),'r')
  hold off
end
title( max(A) );
drawnow
