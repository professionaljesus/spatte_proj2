function [C,G,G2,FV]=matern_prec_matrices(U,FV)
% MATERN_PREC_MATRICES Calculate matrices need to build Matérn precisions
%
%  [C,G,G2]=matern_prec_matrices(U)
% or
%  [C,G,G2]=matern_prec_matrices(U,FV)
%
%  U: N-by-d matrix of d-dimensional coordinate row vectors.
%     Column 1 holds the 1st dimension coordinates, etc.
%  C,G,G2: sparse N-by-N matrices, that can be used to construct
%          precision matrices for Matérn/SPDE fields.
%  FV: A triangular graph specification (optional on R^2, mandatory in R^3)
%
%  Example, alpha=2, range=10:
%   sz = [60,80];
%   [u1,u2] = ndgrid(1:sz(1),1:sz(2));
%   [C,G,G2] = matern_prec_matrices([u1(:),u2(:)]);
%   kappa = 0.3;
%   Q = kappa^4*C + 2*kappa^2*G + G2;

% $Id: matern_prec_matrices.m 4836 2014-12-10 11:09:32Z johanl $

if (nargin<2), FV = []; end

N = size(U,1);
d = size(U,2);

if (d==1)
  error(sprintf('d=%d has not been implemented yet.',d))
end

issphere = 0;
if isempty(FV)
  if (d>2)
    error('For d>2, FV must be supplied.')
  else
    % Avoid colinearities on the boundary:
    U1_ = (max(U(:,1))+min(U(:,1)))/2;
    U2_ = (max(U(:,2))+min(U(:,2)))/2;
    U1 = (U(:,1)-min(U(:,1)))/(max(U(:,1))-min(U(:,1)));
    U2 = (U(:,2)-min(U(:,2)))/(max(U(:,2))-min(U(:,2)));
    U1__ = U1+(U1-0.5).*sin(pi*U2)*1e-6;
    U2__ = U2+(U2-0.5).*sin(pi*U1)*1e-6;
    FV = delaunay(U1__,U2__);
    %    U = [U1__*(max(U(:,1))-min(U(:,1)))+min(U(:,1)),...
    %         U2__*(max(U(:,2))-min(U(:,2)))+min(U(:,2))];
    %    plot(U(:,1),U(:,2),'.')
    %    hold on
    %    trisurf(FV,U(:,1),U(:,2),U(:,1)*0)
    %    hold off
  end
else
  if (d==3)
    issphere = (norm(sum(U.^2,2)-1,inf)<1000*eps);
  end
end

[C,G,G2] = hilbert_appr_sphere(FV,U',issphere);


function [C,G,G2,Cexact]=hilbert_appr_sphere(FV,P,sphere)

nV = size(P,2);
nF = size(FV,1);

Ii = zeros(nF*3,3);
Ij = zeros(nF*3,3);

Cz = zeros(nV,1);
Gz = zeros(nF*3,3);

for f=1:nF
  r = P(:,FV(f,[3,1,2]))-P(:,FV(f,[2,3,1]));
  if (sphere)
    P_ = P(:,FV(f,[1,2,3]));
    s = P_'*P_;
    theta(1) = acos((s(2,3)-s(1,2)*s(1,3))/sqrt((1-s(1,2)^2)*(1-s(1,3)^2)));
    theta(2) = acos((s(3,1)-s(2,3)*s(2,1))/sqrt((1-s(2,3)^2)*(1-s(2,1)^2)));
    theta(3) = acos((s(1,2)-s(3,1)*s(3,2))/sqrt((1-s(3,1)^2)*(1-s(3,2)^2)));
    f_area = sum(theta)-pi; 
  else 
    r1 = r(:,1); r2 = r(:,2);
    f_area = sqrt((r1'*r1)*(r2'*r2)-(r1'*r2)^2)/2; %area of triangle
    if (abs(imag(f_area))>0)
      warning('complex area!')
    end
  end
  Cz(FV(f,:)) = Cz(FV(f,:))+f_area/3;
  Cexactz(3*(f-1)+(1:3),:) = f_area/12*[2 1 1;1 2 1;1 1 2];
  
  Ii(3*(f-1)+(1:3),:) = FV(f,:)'*ones(1,3);
  Ij(3*(f-1)+(1:3),:) = ones(3,1)*FV(f,:);
  Gz(3*(f-1)+(1:3),:) = (r'*r)/(4*f_area);
end

C = spdiags(Cz,0,nV,nV);
G = sparse(Ii(:),Ij(:),Gz(:),nV,nV);
G2 = spdiags(1./sqrt(Cz),0,nV,nV)*G;
G2 = G2'*G2;
Cexact = sparse(Ii(:),Ij(:),Cexactz(:),nV,nV);
