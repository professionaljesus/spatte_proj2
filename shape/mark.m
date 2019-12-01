function px=mark(x,N)
% MARK Mark landmarks in an image.
%
%  marks=mark(x)
%  marks=mark(x,p)
%
%  Mark points with the left mouse button.
%  Press any other mouse button to abort.
%
%  Warning: Since 2004-11-29, a p-by-2 matrix is returned,
%                             and not a 2-by-p matrix.
%           The old behaviour is obtained by using the version in
%           deformable/mark.m instead of deform/mark.m
%           Sorry for the inconvenience.

% Copyright (c) 2001-2002,2004 by Finn Lindgren
% $Revision: 2957 $  $Date: 2006-09-25 09:20:23 +0200 (mån, 25 sep 2006) $

if nargin<2
  N=inf;
end

if size(x,3)>1
  image(x)
  axis xy
else
  imagesc(x)
  axis xy
end

if N<Inf
  px=zeros(N,2);
else
  px=zeros(0,2);
end
for k=1:N
  [u,v,btn] = ginput(1);
  if (btn==1)
    px(k,2) = v;
    px(k,1) = u;
    hold on
    h=text(px(k,1),px(k,2),int2str(k));
    set(h,'color','c','horizontal','center','vertical','middle')
%    plot(px(k,1),px(k,2),'oc')
    hold off
  else
    px(k:end,:) = [];
    return
  end
end
