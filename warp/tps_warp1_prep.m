function prep1=tps_warp1_prep(varargin)
% TSP_WARP1_PREP  Precompute data for use in  tps_warp1
%
%  prep1=tsp_warp1_prep(p0,p1,sz0)
%  prep1=tsp_warp1_prep(p0,p1,sz0,sz1)
%  prep1=tsp_warp1_prep(p0,p1,sz0,sz1,axis0,axis1)
%  prep1=tsp_warp1_prep(p0,p1,sz0,sz1,axis0,axis1,lambda)
%  prep1=tsp_warp1_prep(...,plotflag)
%
%  p0 = The zero-energy control point configuration, np-by-2 matrix.
%  p1 = The target control point configuration, np-by-2 matrix.
%  sz0 = [ny0 nx0] = number of rows, number of columns in X0
%  sz1 = [ny1 nx1] = number of rows, number of columns in X1
%                    Default: sz1 = sz0
%  axis0, axis1: [xmin xmax ymin ymax] for the image axes
%                Default: axis0 = [1 nx0 1 ny0]
%                         axis1 = [1 nx1 1 ny1]
%  lambda: smoothing parameter.
%          lambda=0 (default) gives an interpolating spline.
%  plotflag: If plotflag=1, plot information during the numerical
%            inversion process.
%            Default: plotflag = 0;

% $Id: tps_warp1_prep.m 2957 2006-09-25 07:20:23Z johanl $

p0 = varargin{1};
p1 = varargin{2};
sz0 = varargin{3};
if (nargin<4), sz1 = sz0; else, sz1 = varargin{4}; end
if (nargin<5), axis0 = []; else, axis0 = varargin{5}; end
if (nargin<6), axis1 = []; else, axis1 = varargin{6}; end
if (nargin<7), lambda = []; else, lambda = varargin{7}; end
if (nargin<8), plotflag = []; else, plotflag = varargin{8}; end

if isempty(axis0), axis0 = [1 sz0(2) 1 sz0(1)]; end
if isempty(axis1), axis1 = [1 sz1(2) 1 sz1(1)]; end
if isempty(lambda), lambda = 0; end
if isempty(plotflag), plotflag = 0; end

ny0 = sz0(1);
nx0 = sz0(2);
n0 = ny0*nx0;
ny1 = sz1(1);
nx1 = sz1(2);
n1 = ny1*nx1;
np = size(p0,1);

[x0,y0] = meshgrid(linspace(axis0(1),axis0(2),nx0),...
                   linspace(axis0(3),axis0(4),ny0));
prep0 = tps_warp0_prep(tps_warp0_prep(tps_warp0_prep(...
    p0,sz0,axis0,lambda),x0,y0),p1,sz1,axis1);

if plotflag
  figure(1)
  subplot(221)
  mesh(prep0.x0,prep0.y0,zeros(ny0,nx0))
  view(0,90),axis xy,axis equal
  ax0 = axis;
  subplot(222)
  tmp = reshape(prep0.B*p1,[ny0,nx0,2]);
  mesh(tmp(:,:,1),tmp(:,:,2),zeros(ny0,nx0))
  view(0,90),axis xy,axis equal
  ax1 = axis;
end

%[x1,y1] = meshgrid(1:nx1,1:ny1);
[x1,y1] = meshgrid(linspace(axis1(1),axis1(2),sz1(2)),...
                   linspace(axis1(3),axis1(4),sz1(1)));
xy1 = [x1(:) y1(:)];

xy0 = xy1;

p1_ = p1;
loop = 1;
xy0_(:,:,loop) = xy0;
%jj = [0,logspace(-1,0,8)];
%jj = [0,0.25,0.5,0.75,1];
jj = [0,1];
for j=2:length(jj)
  loop = loop+1;
  p1 = p0*(1-jj(j))+jj(j)*p1_;

% Calculate warp coefficients:
wa = prep0.L\[p1;zeros(3,2)];
w = wa(1:end-3,:);
a = wa(end-2:end,:);

% Initial estimate:
if (j==2)
  %% Use no information:
%  xy0 = xy1;
  %% Use only the affine transformation:
  %% [x1;y1] = a'*[1;x0,y0]
  %% [x1,y1]-a(1,:) = [x0,y0]*a(2:3,:)
%  xy0 = (xy1-a(ones(n,1),:))/a(2:3,:);
  %% Use the reverse TPS as approximation:
  prep_rev = tps_warp0_prep(tps_warp0_prep(p1,sz1,axis1,prep0.lambda),...
                            xy1(:,1),xy1(:,2));
  xy0 = prep_rev.B*p0;
elseif (j==3)
  xy0 = xy0_(:,:,loop-1)+(xy0_(:,:,loop-1)-xy0_(:,:,loop-2))/...
	(jj(loop-1)-jj(loop-2))*(jj(loop)-jj(loop-1));
elseif (j>3)
  delta2 = (xy0_(:,:,loop-1)-xy0_(:,:,loop-2))/(jj(loop-1)-jj(loop-2));
  delta3 = (xy0_(:,:,loop-1)-xy0_(:,:,loop-3))/(jj(loop-1)-jj(loop-3));
  b = (delta3-delta2)/(jj(loop-3)-jj(loop-2));
  a = delta2 + (jj(loop-1)-jj(loop-2))*b;
  xy0 = xy0_(:,:,loop-1)+(jj(loop)-jj(loop-1))*a+...
	(jj(loop)-jj(loop-1))^2*b;
end

[xy1_residual,xy1_grad]=resid(xy0,xy1,p0,wa);

if (plotflag)
  subplot(223)
  mesh(reshape(xy0(:,1),[ny1,nx1]),reshape(xy0(:,2),[ny1,nx1]),zeros(ny1,nx1))
  view(0,90),axis xy,axis equal,axis(ax0)
  title('Initial guess')
  subplot(224)
  mesh(reshape(xy1(:,1)+xy1_residual(1:n1,1),[ny1,nx1]),...
       reshape(xy1(:,2)+xy1_residual(n1+1:end,1),[ny1,nx1]),zeros(ny1,nx1))
  view(0,90),axis xy,axis equal,axis(ax1)
  title([sqrt(xy1_residual'*xy1_residual),...
	 max(sqrt(xy1_residual(1:n1).^2+xy1_residual(n1+1:end).^2))])
  drawnow
  %pause
end

% Optimise:
innerloop = 0;
while ((innerloop<10) &...
       (max(xy1_residual(1:n1).^2+xy1_residual(n1+1:end).^2) >= 0.01.^2))
  innerloop = innerloop+1;
%  figure(2)
  updated = logical(zeros(n1,1));
  for i=1:n1
    if ((xy1_residual(i).^2+xy1_residual(n1+i).^2) >= 0.01.^2)
      updated(i) = logical(1);
      % zeros(2,1) = xy1_residual(i+[0,n1]) + xy1_grad(i+[0,n1],:) * h
      h = - xy1_grad(i+[0,n1],:) \ xy1_residual(i+[0,n1]);
      len = sqrt(h'*h);
%      if (len>10)
%	h = h/len*5;
%      end
      len = sqrt(h'*h);
      Q(1:2,1) = xy1_residual(i+[0,n1]);

      k = 0; step = 1; steps = 0;
      finished = 0;
      while (~finished & (k<10))
	k = k+1;
	steps(1+k) = step;
	xy0_new = xy0(i,:) + h'*step;
	Q(1:2,1+k) = resid(xy0_new,xy1(i,:),p0,wa);
	finished = (sum(Q(:,1+k).^2)<sum(Q(:,1).^2));
	step = step*0.5;
      end
      if (k>50) & plotflag
	figure(2)
	[steps_,idx_] = sort(steps);
	plot(steps_,sqrt(sum(Q(:,idx_).^2,1)),'r',...
	     steps_,abs(Q(1,idx_)),'b',...
	     steps_,abs(Q(2,idx_)),'g')
	title([k,step*2])
	drawnow
	figure(1)
%	pause
      end
      xy0(i,:) = xy0_new;
    end
  end
%  pause
%  figure(1)
  

  [xy1_residual_,xy1_grad_]=resid(xy0(updated,:),xy1(updated,:),p0,wa);
  xy1_residual([updated;updated],:) = xy1_residual_;
  xy1_grad([updated;updated],:) = xy1_grad_;


  if (plotflag)
    subplot(222)
    tmp = reshape(prep0.B*p1,[ny0,nx0,2]);
    mesh(tmp(:,:,1),tmp(:,:,2),zeros(ny0,nx0))
    view(0,90),axis xy,axis equal,axis(ax1)
    subplot(223)
    mesh(reshape(xy0(:,1),[ny1,nx1]),reshape(xy0(:,2),[ny1,nx1]),zeros(ny1,nx1))
    view(0,90),axis xy,axis equal,axis(ax0)
    title(['Iteration ' int2str(innerloop)])
    subplot(224)
    mesh(reshape(xy1(:,1)+xy1_residual(1:n1,1),[ny1,nx1]),...
	 reshape(xy1(:,2)+xy1_residual(n1+1:end,1),[ny1,nx1]),zeros(ny1,nx1))
    view(0,90),axis xy,axis equal,axis(ax1)
    title([sqrt(xy1_residual'*xy1_residual),...
	   max(sqrt(xy1_residual(1:n1).^2+xy1_residual(n1+1:end).^2))])
    drawnow
  end

%pause

end

xy0_(:,:,loop) = xy0;

if plotflag
  disp(xy1_residual'*xy1_residual)

  figure(2)
  clf
  plot(jj(1:j),reshape(xy0_(ny1*floor(nx1/8)+(1:ny1),1,:),[ny1,j]))
  figure(1)
end

end

if plotflag
  for loop=2:length(jj)
    [norm((xy0_(:,1,loop)-xy1(:,1)) - (xy0(:,1)-xy1(:,1))),...
     norm((xy0_(:,1,loop)-xy1(:,1))/jj(loop) - (xy0(:,1)-xy1(:,1)))]
  end
end

prep1.prep0 = prep0;
prep1.xy0 = xy0;


function [xy1_residual,xy1_grad]=resid(xy0,xy1,p0,wa)

np = size(p0,1);
n = size(xy0,1);
  
x0d = (xy0(:,ones(1,np))-ones(n,1)*p0(:,1)');
y0d = (xy0(:,ones(1,np)+1)-ones(n,1)*p0(:,2)');
r2 = x0d.^2+y0d.^2;
lr2 = log(r2+(r2==0));
U = (r2.*lr2)/2;

tmp = [U,[ones(n,1),xy0]];
xy1_residual_x = tmp*wa(:,1) - xy1(:,1);
xy1_residual_y = tmp*wa(:,2) - xy1(:,2);
xy1_residual = [xy1_residual_x;xy1_residual_y];

if (nargout>1)
  lr2 = max(1,lr2);
%  lr2 = lr2*1;
  % x*(1+ln(x^2+y^2)) = x*()

  xy1_grad_x = [[x0d.*(1+lr2),...
		 zeros(n,1),ones(n,1),zeros(n,1)]*wa(:,1),...
		[y0d.*(1+lr2),...
		 zeros(n,1),zeros(n,1),ones(n,1)]*wa(:,1)];
  xy1_grad_y = [[x0d.*(1+lr2),...
		 zeros(n,1),ones(n,1),zeros(n,1)]*wa(:,2),...
		[y0d.*(1+lr2),...
		 zeros(n,1),zeros(n,1),ones(n,1)]*wa(:,2)];
  xy1_grad     = [xy1_grad_x;xy1_grad_y];
end
