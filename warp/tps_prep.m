function prep=tps_prep(varargin)
% TPS_PREP  Precompute data for use in  tps_pull  and tps_push
%
%  prep=tps_prep(p0,lambda)
%  prep=tps_prep(prep,p1)
%  prep=tps_prep(prep,sz0,axis0)
%  prep=tps_prep(prep,xy0)
%
%  In the second form,  prep  is a struct given by a previous
%  invocation of  tps_prepare.  Both calls must be used before
%  using prep as input to tps_pull or tps_push.
%
%  p0 = The zero-energy control point configuration, np-by-2 matrix.
%  lambda = smoothing parameter.
%           lambda=0 gives an interpolating spline.

% $Id: tps_prep.m 3176 2006-11-27 07:49:23Z johanl $

%% Determine the call type:
if (~isstruct(varargin{1}))
  call_type = 1;
elseif (nargin>2)
  call_type = 3;
elseif (nargin<2)
  error('Too few input arguments.')
elseif (size(varargin{2},1)==1)
  call_type = 3;
elseif (size(varargin{2},3)>1)
  call_type = 4;
else
  call_type = 2;
end

switch (call_type)
 case 1,
  p0 = varargin{1};
  if (nargin<2), lambda = []; else lambda = varargin{2}; end
  if isempty(lambda), lambda = 0; end
  prep.p0 = p0;
  prep.lambda = lambda;
  
  np = size(p0,1);
  prep.np = np;
  
  p0x = p0(:,ones(1,np));
  p0y = p0(:,1+ones(1,np));
  r2 = (p0x-p0x').^2 + (p0y-p0y').^2;
  K = r2.*log(r2+(r2==0))/2;
  P = [ones(np,1), p0];
  prep.L = [K+eye(np)*lambda, P; P', zeros(3,3)];
  
  prep.available_prep = [0,0];

 case 2,
  prep = varargin{1};
  prep.p1 = varargin{2};
  prep.available_prep(1) = 1;
  
  prep.Lp1 = prep.L\[prep.p1;zeros(3,2)];

  if (prep.available_prep(2))
    prep.xy1 = reshape(prep.Upoints*prep.Lp1,[prep.sz0+1,2]);
  end

 case 3,
  prep = varargin{1};
  sz0 = varargin{2};
  if (nargin<3), axis0 = []; else axis0 = varargin{3}; end
  if isempty(axis0), axis0 = [1 sz0(2) 1 sz0(1)]+[-0.5,0.5,-0.5,0.5]; end
  prep.sz0 = sz0;
  prep.axis0 = axis0;
  [x0,y0] = ...
      meshgrid(linspace(axis0(1),axis0(2),sz0(2)+1),...
               linspace(axis0(3),axis0(4),sz0(1)+1));
  prep.xy0 = x0;
  prep.xy0(:,:,2) = y0;
 
 case 4,
  prep = varargin{1};
  prep.xy0 = varargin{2};
  prep.sz0 = [size(prep.xy0,1),size(prep.xy0,2)]-1;
  x0 = prep.xy0(:,:,1);
  y0 = prep.xy0(:,:,2);

end

if ((call_type==3) || (call_type==4))
  x0 = x0(:);
  y0 = y0(:);
  n = (prep.sz0(1)+1)*(prep.sz0(2)+1);
  p0 = prep.p0';
  r2 = (x0(:,ones(1,prep.np))-ones(n,1)*p0(1,:)).^2 +...
       (y0(:,ones(1,prep.np))-ones(n,1)*p0(2,:)).^2;
  U = r2.*log(r2+(r2==0))/2;
  prep.Upoints = [U,ones(n,1),x0,y0];

  prep.available_prep(2) = 1;

  if (prep.available_prep(1))
    prep.xy1 = reshape(prep.Upoints*prep.Lp1,[prep.sz0+1,2]);
  end

end
