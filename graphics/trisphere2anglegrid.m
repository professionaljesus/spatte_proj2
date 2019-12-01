function [gridval,FF]=trisphere2anglegrid(FV,S,val,u,v,FF)
% TRISPHERE2ANGLEGRID Map a spherical field to a flat projection
%
% gridval=trisphere2anglegrid(FV,S,val,u,v)
% FV contains the triangulation
% S cointains the 3D-points
% val are the a values at each point in S 
%     (or for each triangle in FV).
% u,v is a longitude & latitude grid in radians defining at what
%     resolution to compute the flat projection
%
% Example, 2-degree spacing:
%  u = linspace(-180,180,181);
%  v = linspace(-90,90,91);
%  gridval=trisphere2anglegrid(FV,S,val,u*pi/180,v*pi/180);
% Plot with surf:
%  surf(u,v,gridval*0,gridval)
%  shading interp,axis([-180,180,-90,90])
% Plot with imagesc:
%  imagesc(u,v,gridval)
%  axis xy

% $Id: trisphere2anglegrid.m 4586 2012-10-08 16:18:33Z johanl $

if (nargin<6), FF = []; end
%check orientation of FV
v1 = S(FV(:,2),:) - S(FV(:,1),:);
v2 = S(FV(:,3),:) - S(FV(:,2),:);
tmp = sum(cross(v1,v2) .* S(FV(:,1),:),2);
%reorder triangles that are oriented "inwards"
FV(tmp<0,:) = FV(tmp<0,[3 2 1]);

if isempty(FF), FF = fv2ff(FV); end

if (min(size(u))==1) && (min(size(v))==1)
  [u,v] = meshgrid(u,v);
elseif any(size(u)~=size(v))
  error('u and v must be vectors, or have the same size')
end
% Transpose u and v for more efficient triangle location:
u = u'; v = v';

if (size(val,1)==size(FV,1))
  Fval = true;
  valdim = size(val,2);
elseif (size(val,1)==size(S,1))
  Fval = false;
  valdim = size(val,2);
else
  error('The size of val must match the length of FV or S')
end
gridval = zeros([numel(u),valdim]);

f = 1;
for k=1:length(u(:))
  P0 = [[cos(u(k));...
         sin(u(k))]*cos(v(k));...
        sin(v(k))];
  [f,w] = point_locate(FV,S',P0,f,1,FF,0);
  if (Fval)  
    gridval(k,:) = val(f,:);
  else
    gridval(k,:) = w'*val(FV(f,:),:);
  end
end

gridval = reshape(gridval,[size(u),valdim]);

u = u';
v = v';
gridval = permute(gridval,[2,1,3]);

if (0)
surf(u,v,gridval(:,:,1))
shading interp
view(0,-90)
axis tight
drawnow
end



function [f,w,f_seq_len,f_sequence]=point_locate(FV,P,pt,f,is_sphere,FF,plotflag)

if (nargin<4), f = []; end
if (nargin<5), is_sphere = []; end
if (nargin<6), FF = []; end
if (nargin<7), plotflag = []; end
if isempty(f), f = ceil(rand(1,1)*size(FV,1)); end
if isempty(is_sphere), is_sphere = false; end
if isempty(FF), FF = fv2ff(FV); end
if isempty(plotflag), plotflag = false; end

d = size(P,1);
if (is_sphere && (d~=3))
  error('A sphere must be embedded in R^3, not R^%i.',d)
end

p_search = pt;
f_sequence = f;

not_found = true;
while (not_found)
  PT = P(:,FV(f,:));
  ET = PT(:,[3,1,2])-PT(:,[2,3,1]);
  if (d==3)
    n = cross_int(ET(:,1),ET(:,2));
    n = n/sqrt(n'*n);
  end
  if (is_sphere)
    %    if (mean(PT,2)'*n<0)
    %   n = -n;
    %   warning('Flipped normal!')
    % end
    cos_theta = pt'*n;
    if (cos_theta <= 0)
      p_search = pt-cos_theta*n;
      p_search = p_search/sqrt(p_search'*p_search)+n;
    else
      p_search = pt;
    end
    w = [PT;ones(1,3)]\[p_search*(PT(:,1)'*n)/(p_search'*n);1];
  else
    if (d==2)
      w = [PT;ones(1,3)]\[p_search;1];
    else
      w = [PT,n;ones(1,3),0]\[p_search;1];
      w = w(1:3);
    end
  end
  % Th = [ET(:,2)-ET(:,1)*(ET(:,2)'*ET(:,1))/(ET(:,1)'*ET(:,1)),...
  %       ET(:,3)-ET(:,2)*(ET(:,3)'*ET(:,2))/(ET(:,2)'*ET(:,2)),...
  %       ET(:,1)-ET(:,3)*(ET(:,1)'*ET(:,3))/(ET(:,3)'*ET(:,3))];
  %  distance = w.*sqrt(sum(Th.^2,1))';
  [ws,i] = min(w);
  if (w(i(1)) >= -eps*10)
    not_found = false;
  else
    f_ = FF(f,rem(i(1),3)+1);
    if (any(f_sequence==f_))
      not_found = false;
      if (min(w)<-eps*100)
        warning('It seems we got stuck somewhere! w_min=%f',min(w)/eps)
      end
    else
      f = f_;
      f_sequence = [f_sequence,f];
    end
  end
  if plotflag && (~not_found)
    %    trimesh(FV,P(1,:),P(2,:),P(3,:))
    trimesh(FV(f_sequence,:),P(1,:)*1.01,P(2,:)*1.01,P(3,:)*1.01,...
            'edgecolor','r')
    hold on
    PT_ = mean(PT,2);
    plot3(PT_(1),PT_(2),PT_(3),'g*')
    plot3(pt(1),pt(2),pt(3),'b*')
    plot3(p_search(1),p_search(2),p_search(3),'r*')
    hold off
    axis([-1 1 -1 1 -1 1]*1.5)
    view(0,90)
    drawnow
    pause(0)
  end
end

f_seq_len = length(f_sequence);

function C=cross_int(A,B)
C = [A(2)*B(3) - A(3)*B(2);...
     A(3)*B(1) - A(1)*B(3);...
     A(1)*B(2) - A(2)*B(1)];



function FF=fv2ff(FV,method)
% FV2FF Calculates FF graph from FV graph.
%
%  CALL: FF=fv2ff(FV)

% History since 1998-05-27
% 1999-06-29: Added 'uint'-tolerance to method 1

if nargin<2
  if (exist('ismembc','file')>0)
    method = 1;
  else
    method = 2;
  end
end

% New algorithm Nr2 1998-05-26
% True linear nF complexity.
% Faster than the 05-25 algorithm for all nF.
% Example: nP=10000, nF=19973
% Example: nP=10000, nF=19973
%          T_new2   =   17 sec (with direct call to ismembc, 45 otherwise)
%          T_new    =   54 sec
%          T_old    =  913 sec
%          T_oldest ~ 2700 sec
%          (Relative speeds: 1:3:50:150)
switch method
 case 1,
  FV=double(FV);
  nF=size(FV,1);
  nP=max(FV(:));
  VFs=sparse(FV(:),[1:nF 1:nF 1:nF]',1,nP,nF);
  FFs=(VFs'*VFs==2);
  [i,j]=find(FFs);
  FVsort=sort(FV,2); % For use with ismembc
  FF=zeros(nF,3);
  for k=1:length(i)
    ik=i(k);
    jk=j(k);
    ii=find(ismembc(FV(ik,:),FVsort(jk,:)));
    
    if ii(1)+1==ii(2)
      ii=ii(1);
    else
      ii=ii(2);
    end
    FF(ik,ii)=jk;
    FF(jk,FV(jk,:)==FV(ik,rem(ii,3)+1))=ik;
  end
  
   
   % New algorithm 1998-05-26
   % Approx nF^1--2 complexity.
   % Faster than the 05-25 algorithm for all nF.
 case 2,
  nF=size(FV,1);
  nP=max(FV(:));
  VFs=sparse(FV(:),[1:nF 1:nF 1:nF]',1,nP,nF);
  FFs=(VFs'*VFs==2);
  [i,j]=find(FFs);
  FF=zeros(nF,3);
  for k=1:length(i)
    ik=i(k);
    jk=j(k);
    % Important: The complexity of find for sparse matrices
    % is much lower for columns than for rows.
    v=find(VFs(:,ik).*VFs(:,jk))';
    ii=1;
    if (FV(ik,ii)==v(1)) && (FV(ik,ii)==v(2))
      ii=2;
      if (FV(ik,ii)==v(1)) && (FV(ik,ii)==v(2))
	ii=3;
      end
    end
    ii=rem(ii,3)+1;
    jj=1;
    if (FV(jk,jj)==v(1)) && (FV(jk,jj)==v(2))
      jj=2;
      if (FV(jk,jj)==v(1)) && (FV(jk,jj)==v(2))
	jj=3;
      end
    end
    jj=rem(jj,3)+1;
    FF(ik,ii)=jk;
    FF(jk,jj)=ik;
  end
  
  
  
  %% New algorithm 1998-05-25
  %% Relative speed new:old approx 1:3
  %% Approx nF^2 complexity.
 case 3
  nF=size(FV,1);
  FF=zeros(nF,3);
  for i=1:3
    j=rem(i,3)+1;
    for f=1:nF
      p=FV(f,[1 2 3 1]);
      FVi=FV(f+1:end,i);
      FVj=FV(f+1:end,j);
      for k=1:3
	f2=find(((FVj==p(k)).*(FVi==p(k+1))))+f;
	if ~isempty(f2)
	  FF(f,k)=f2;
	  FF(f2,i)=f;
	end
      end
    end
  end
  
  % Old algorithm:
  % Approx >nF^2 complexity.
 case 4,
  nF=size(FV,1);
  FF=zeros(nF,3);
  for f=1:nF
    p=FV(f,[1 2 3 1]);
    for k=1:3
      f2=setdiff(...
	  find((FV(:,2)==p(k) & FV(:,1)==p(k+1)) |...
	       (FV(:,3)==p(k) & FV(:,2)==p(k+1)) |...
	       (FV(:,1)==p(k) & FV(:,3)==p(k+1))),...
	  f);
      if ~isempty(f2)
	FF(f,k)=f2;
      end
    end
  end
end
