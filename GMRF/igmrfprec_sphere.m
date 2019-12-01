function [Q,W,area,W1,pole_weight]=igmrfprec_sphere(sz,direction)
% IGMRFPREC_SPHERE Constructs a precision matrix for spherical IGMRF
%
%   Q = igmrfprec_sphere(sz)
%   [Q,W] = igmrfprec_sphere(sz)
%   [Q,W,area] = igmrfprec_sphere(sz)
%
% The intrinsic random field X \in N(0, Q^-1) is an approximate solution
% to the stochastic partial differential equation
%   Laplace(x(u)) = spatial white noise with variance 1
% on the unit radius sphere.
% The generated model will only give approximately equal values for
% the grid points all corresponding to the north pole (and the south
% pole), but the practical differences are small.
%
%   sz : The grid size,  sz = [m,n] for an m-by-n image
%        sz(1) is the number of latitudes
%              row 1 = north pole = +90 degrees,
%              row m = south pole = -90 degrees
%        sz(2) is the number of longitudes
%              column 1 = -180 degrees,
%              column n = +180-360/sz(2) degrees
%              column n/2 = 0 degrees, the Greenwich meridian.
%
%   Q    : The sparse precision matrix
%   W    : The Laplacian increment matrix: Q=W'*diag(area)*W
%          Unlike for the igmrfprec function, the increment definitions
%          are different for every grid row, to compensate for the
%          difference in distances on the sphere for different latitudes.
%   area : The local area around each grid node, stored in colstack mode.

% $Id: igmrfprec_sphere.m 4836 2014-12-10 11:09:32Z johanl $

if (nargin<2), direction = []; end
if (isempty(direction)), direction = 0; end

pole_weight = 1e6;

lat = linspace(90,-90,sz(1))/360*2*pi;
lon = linspace(-180,180-360/sz(2),sz(2))/360*2*pi;

P1 = cos(lat)'*cos(lon);
P2 = cos(lat)'*sin(lon);
P3 = sin(lat)'*ones(1,sz(2));

P = [P1(:),P2(:),P3(:)];

% (I,J) (I+1,J) (I+1,J+1)  K00 K10 K11
% (I,J) (I+1,J+1) (I,J+1)  K00 K11 K01
K00 = (1:sz(1)-1)'*ones(1,sz(2))+...
      sz(1)*ones(sz(1)-1,1)*(0:sz(2)-1);
K10 = K00+1;
K01 = rem(K00+sz(1)-1,sz(1)*sz(2))+1;
K11 = K01+1;
if (direction==0)
  FV = [K00(:) K10(:) K11(:);...
        K00(:) K11(:) K01(:)];
else
  FV = [K00(:) K10(:) K01(:);...
        K10(:) K11(:) K01(:)];
end
for k=1+sz(1)*(1:sz(2)-1)
  FV(FV==k) = 1;
  FV(FV==sz(1)*sz(2)-k+1) = sz(1)*sz(2);
end
FV((FV(:,1)==FV(:,2))|(FV(:,2)==FV(:,3))|(FV(:,1)==FV(:,3)),:) = [];

[W,Av] = laplacian_sphere_weights(FV,P',[]);

area = Av';

II1 = [];
JJ1 = [];
KK1 = [];

i = 1;
for j=1:sz(2)
  II1 = [II1;...
         [1;1]*(i+sz(1)*(j-1))];
  JJ1 = [JJ1;...
         [i+sz(1)*(j-1);i+sz(1)*(mod(j,sz(2)))]];
  KK1 = [KK1;...
         [1;-1]];
end

i = sz(1);
for j=1:sz(2)
  II1 = [II1;...
         [1;1]*(i+sz(1)*(j-1))];
  JJ1 = [JJ1;...
         [i+sz(1)*(j-1);i+sz(1)*(mod(j,sz(2)))]];
  KK1 = [KK1;...
         [1;-1]];
end

W1 = sparse(II1,JJ1,KK1,sz(1)*sz(2),sz(1)*sz(2));

Q = W1'*W1*pole_weight + W'*spdiags(area,0,length(W),length(W))*W;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,Av]=laplacian_sphere_weights(FV,P,filenameroot,method)
% laplacian_sphere_weights(FV,P,'filenameroot',method)
% produces files filenameroot.graph and filenameroot.weights
% for input to GMRFLib.
% Assumes that VV,P describes a triangulated sphere.
%
% Example:
%   [FV,P]=trisphere(5);
%   laplacian_sphere_weights(FV,P,'sphere5');
%
% Uses laplacian_weights_calc.m

% Finn Lindgren 2004
% $Id: igmrfprec_sphere.m 4836 2014-12-10 11:09:32Z johanl $

if (nargin<4), method = []; end
if (nargin<3), filenameroot = []; end

if isempty(method), method = 'poly'; end

At = zeros(size(FV,1),1);
Av = zeros(1,size(P,2));
for t=1:size(FV,1)
  At(t) = spherical_triangle_area(P(:,FV(t,:)));
  Av(FV(t,:)) = Av(FV(t,:)) + At(t)/3;
end

VV = fv2vv(FV);

switch (method)
 case 'poly',
  W = calc_poly(FV,VV,P,At,Av);
 case 'sh',
  W = calc_sh(FV,VV,P,At,Av);
 case 'fa',
  W = calc_fa(FV,VV,P,At,Av);
 case 'fa2',
  W = calc_fa2(FV,VV,P,At,Av);
 otherwise,
  error(['Unknown method "',method,'".']);
end

if (~isempty(filenameroot))
  print_weights(FV,VV,P,W,Av,filenameroot);
end


%%%%%%%%
function W=calc_poly(FV,VV,P,At,Av)

W = sparse(size(P,2),size(P,2));
for v=1:size(P,2)
  Nv = find(VV(v,:));
  d = length(Nv);
  if (d==0), continue, end % No neighbours.
  % Make an (arbitrary) orthogonal coordinate system, with P(:,v) as the
  % first vector:
  B = [P(:,v), rand(3,2)];
  B(:,2) = B(:,2) - P(:,v)*(B(:,2)'*P(:,v));
  B(:,2) = B(:,2)/sqrt(B(:,2)'*B(:,2));
  B(:,3) = B(:,3) - P(:,v)*(B(:,3)'*P(:,v));
  B(:,3) = B(:,3) - B(:,2)*(B(:,3)'*B(:,2));
  B(:,3) = B(:,3)/sqrt(B(:,3)'*B(:,3));
  % Make it right-handed:
  if (det(B)<1), B(:,3) = -B(:,3); end
  % Compute the coordinates:
  Ploc = B\P(:,Nv);
  
  % Rescale to match spherical arc distance:
%  Ploc = Ploc*diag(acos(P(:,v)'*P(:,Nv))./sqrt(sum(Ploc.^2,1)));
%  Ploc = Ploc*diag(sqrt(sum((P(:,Nv)-P(:,v)*...
%                             ones(1,length(Nv))).^2,1)) ./ ...
%                   sqrt(sum(Ploc.^2,1)));

  Ploc = Ploc*diag(acos(P(:,v)'*P(:,Nv)) ./ ...
                   sqrt(sum((P(:,Nv)-P(:,v)*...
                             ones(1,length(Nv))).^2,1)));
  
  % Calculate the stencil for the Laplacian operator:
%  [w,w0] = laplacian_weights_calc(Ploc(2:3,:));
  
  w = laplacian_weights_calc([[0;0],Ploc(2:3,:)]);
  w0 = w(1);
  w = w(2:end);
  
  W(v,v) = w0;
  W(v,Nv) = w(:)';
end
  
 
%%%%%%%%
function W=calc_sh(FV,VV,P,At,Av)

W = sparse(size(P,2),size(P,2));
for v=1:size(P,2)
  Nv = find(VV(v,:));
  d = length(Nv);
  % Make an (arbitrary) orthogonal coordinate system, with P(:,v) as the
  % first vector:
  B = [P(:,v), rand(3,2)];
  B(:,2) = B(:,2) - P(:,v)*(B(:,2)'*P(:,v));
  B(:,2) = B(:,2)/sqrt(B(:,2)'*B(:,2));
  B(:,3) = B(:,3) - P(:,v)*(B(:,3)'*P(:,v));
  B(:,3) = B(:,3) - B(:,2)*(B(:,3)'*B(:,2));
  B(:,3) = B(:,3)/sqrt(B(:,3)'*B(:,3));
  % Make it right-handed:
  if (det(B)<1), B(:,3) = -B(:,3); end
  % Compute the coordinates:
  Ploc = B\P(:,[v,Nv]);
  Ploc = [Ploc(2,:);Ploc(3,:);Ploc(1,:)];

  w = zeros(1+d,1);
  V_acc = eye(1+d);
  for o=0:1
    if (o==0)
      A = spherical_harmonics(o,Ploc);
      b = -o*(o+1)*A(:,1);
    elseif (o==-1)
      A = [];
      b = [];
      for oo=1:1
        A_ = spherical_harmonics(oo,Ploc);
        A(end+1:end+2*oo+1,:) = A_*(-oo*(oo+1))^(0-2);
        b(end+1:end+2*oo+1,1) = A_(:,1)*(-oo*(oo+1))^(1-2);
      end
    else
      A = [];
      b = [];
      for oo=1:3 % "3" should be a variable parameter...
        A_ = spherical_harmonics(oo,Ploc);
        A(end+1:end+2*oo+1,:) = A_/(2*oo+1)^4;
        b(end+1:end+2*oo+1,1) = A_(:,1)*(-oo*(oo+1))/(2*oo+1)^4;
      end
    end
    [U,S,V] = svd(A*V_acc);
    kappa = rank(S);
    U1 = U(:,1:kappa);
    U0 = U(:,kappa+1:end);
    S1 = S(1:kappa,1:kappa);
    V1 = V(:,1:kappa);
    V0 = V(:,kappa+1:end);
    w = w + V_acc*(V1*(S1\(U1'*b)));
    V_acc = V_acc*V0;
    if (norm(U0'*b)>1e-7)
%      fprintf('d=%i, LS order: %i\n',d,o);
      break;
    end
  end
  W(v,[v,Nv]) = w';
end




%%%%%%%%
function W=calc_fa(FV,VV,P,At,Av)

D = acos(max(-1,min(1,P'*P)));

W = sparse(size(P,2),size(P,2));
for t=1:size(FV,1)
  W(FV(t,:),FV(t,[2,3,1])) = W(FV(t,:),FV(t,[2,3,1])) +...
      diag(2/3*At(t)./Av(FV(t,:))'./...
           diag(D(FV(t,:),FV(t,[2,3,1]))).^2);
  W(FV(t,[2,3,1]),FV(t,:)) = W(FV(t,[2,3,1]),FV(t,:)) +...
      diag(2/3*At(t)./Av(FV(t,[2,3,1]))'./...
           diag(D(FV(t,:),FV(t,[2,3,1]))).^2);
end
W = W + spdiags(-sum(W,2),0,size(P,2),size(P,2));


%%%%%%%%
function W=calc_fa2(FV,VV,P,At,Av)

W = sparse(size(P,2),size(P,2));
for t=1:size(FV,1)
  A2 = sum(cross(P(:,FV(t,2))-P(:,FV(t,1)),...
                 P(:,FV(t,3))-P(:,FV(t,1))).^2,1);
  d = P(:,FV(t,[3,1,2]))-P(:,FV(t,[2,3,1]));
  B = (d'*d)/A2;
  W(FV(t,:),FV(t,:)) = W(FV(t,:),FV(t,:)) + B*At(t);
end
W = -W./(Av'*ones(1,size(W,1)));


%%%%%%%%
function area=triangle_area(P)
area = cross((P(:,2)-P(:,1)),(P(:,3)-P(:,1)));
area = sqrt(area'*area)/2;

%%%%%%%%
function area=spherical_triangle_area(P)
s = P'*P;
theta(1) = acos((s(2,3)-s(1,2)*s(1,3))/sqrt((1-s(1,2)^2)*(1-s(1,3)^2)));
theta(2) = acos((s(3,1)-s(2,3)*s(2,1))/sqrt((1-s(2,3)^2)*(1-s(2,1)^2)));
theta(3) = acos((s(1,2)-s(3,1)*s(3,2))/sqrt((1-s(3,1)^2)*(1-s(3,2)^2)));
area = sum(theta)-pi;


%%%%%%%%
function print_weights(FV,VV,P,W,Av,filenameroot)

fidG  = fopen([filenameroot,'.graph'],'w');
fidP  = fopen([filenameroot,'.coord'],'w');
fidWA = fopen([filenameroot,'.weights'],'w');
fidW  = fopen([filenameroot,'.laplacian'],'w');
fidA  = fopen([filenameroot,'.areas'],'w');

% Print the number of vertices
fprintf(fidG,'%i\n',size(P,2));

for v=1:size(P,2)
  Nv = find(VV(v,:));
  d = length(Nv);
  % Print the graph and weight structure:
  % v, d, Nv
  fprintf(fidG,['%i %i',repmat(' %3i',[1,d]),'\n'],...
          v-1,d,Nv-1);
  % P
  fprintf(fidP,['%11.16f %11.16f %11.16f','\n'],...
          P(1,v),P(2,v),P(3,v));
  % (w0, w)*sqrt(A)
  fprintf(fidWA,['%12.16f',repmat(' %11.16f',[1,d]),'\n'],...
          full(W(v,v)*sqrt(Av(v))),...
          full(W(v,Nv)*sqrt(Av(v))));
  % w0, w
  fprintf(fidW,['%12.16f',repmat(' %11.16f',[1,d]),'\n'],...
          full(W(v,v)),...
          full(W(v,Nv)));
  % A
  fprintf(fidA,['%12.16f','\n'],...
          full(Av(v)));
end

fclose(fidG);
fclose(fidP);
fclose(fidWA);
fclose(fidW);
fclose(fidA);



function VV=fv2vv(FV,n)
% FV2VV Calculate sparse VV graph matrix from FV graph.

if (nargin<2)
  n = max(FV(:));
end

VV=sparse([FV(:,1);FV(:,1);FV(:,2)],...
          [FV(:,2);FV(:,3);FV(:,3)],...
          ones(size(FV,1)*3,1),n,n);

VV=((VV+VV')>0); % Handle open triangulations by symmetrisation.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,w0]=laplacian_weights_calc(P)
% LAPLACIAN_WIGHT_CALC Compute weights for approximation of the Laplacian.
%
%   w = laplacian_weights_calc(P)
%   P : 2-by-n
%   w : n-by-1
%
% Obsolete syntax:
%   [w,w0] = laplacian_weights_calc(Pn)
%   Pn: 2-by-d, Pn = P(:,2:d+1), P(:,1) is assumed = [0;0]
%   w : d-by-1, the weights for Pn
%   w0: 1-by-1, the weight for [0;0]
%
% For n>=6, and no duplicate points in P, the solution space has
% dimension n-6.
% For n==5 and n==4 only special configurations yield solutions.
% For n<=3 no solutions exist.
% For degenerated cases, an attempt is made to find an approximate
% solution.
% For underdetermined cases, an attempt is made to make the error terms
% isotropic.

% (c) 2004 by Finn Lindgren 

if (nargout>1) % Old syntax.
  w = simple_calc([[0;0],P]);
  w0 = w(1);
  w = w(2:end);
  return;
end

if (size(P,2)>1)
  w = simple_calc(P);
  %w = isotropic_error_calc(P);
else
  w = 0;
end
return;

pr = 0;

n = size(P,2);
theta = zeros(n,1);
VV = eye(n);
least_squares_appr = 0;
if pr, fprintf('Degrees of freedom left:',n); end
for p=0:5
  if (p==0)
    Ap = ones(1,n);
    bp = 0;
    b = bp;
  else
    if (p==1)
      Ap = P;
      bp = [0;0];
    elseif (p==2)
      Ap = [P(1,:).^2;P(1,:).*P(2,:);P(2,:).^2];
      bp = [2;0;2];
%    elseif (p==3)
%      Ap = [P(1,:).^3;P(1,:).^2.*P(2,:);P(1,:).*P(2,:).^2;P(2,:).^3];
%      bp = [0;0;0;0];
%    elseif (p==4)
%      VV = [[VV,zeros(n,1)];zeros(1,size(VV,2)),1];
%      theta = [theta;0];
%      Ap = [[P(1,:).^4;P(1,:).^3.*P(2,:);P(1,:).^2.*P(2,:).^2;...
%             P(1,:).*P(2,:).^3;P(2,:).^4],...
%           [48;0;8;0;48]];
%      bp = [0;0;0;0;0];
    elseif (p==3)
      VV = [[VV,zeros(n,1)];zeros(1,size(VV,2)),1];
      theta = [theta;0];
      Ap = [P(1,:).^3;P(1,:).^2.*P(2,:);P(1,:).*P(2,:).^2;P(2,:).^3];
      Ap = [[Ap;...
            [P(1,:).^4;P(1,:).^3.*P(2,:);P(1,:).^2.*P(2,:).^2;...
             P(1,:).*P(2,:).^3;P(2,:).^4]],...
           [0;0;0;0; 48;0;8;0;48]];
      bp = [0;0;0;0; 0;0;0;0;0];
    elseif (p==4)
      continue;
    elseif (p==5)
      Ap = [zeros(1,n),1];
      bp = [0];
    end
    b = bp-Ap*theta;
  end
  A = Ap*VV;
  [Up,Sp,Vp] = svd(A);
  [m1,m2] = size(A);
  if (min(m1,m2)==1)
    dSp = Sp(1,1);
  else
    dSp = diag(Sp);
  end
  tol = m2*norm(dSp)*eps; % Default tolerance from "rank".
  kp = sum(dSp>tol);
  Ub = Up'*b;
  theta = theta + VV*( Vp(:,1:kp) * (Ub(1:kp)./dSp(1:kp)) );
  if (any(abs(Ub(kp+1:m1))>tol)) % Least squares aproximation.
    least_squares_appr = 1;
  end
  % m2-kp degrees of freedom left.
  if pr, fprintf(' %2i',m2-kp); end
  if (kp>=m2) % No more degrees of freedom left.
    break;
  end
  VV = VV*Vp(:,kp+1:m2);
end
w = theta(1:n);
if (least_squares_appr)
  if pr, fprintf('\tError-order == %i',min(4,p)+(~least_squares_appr)); end
  if pr, fprintf('\tLS-(p,residual): (%i,%f)',p,norm(Ap*theta-bp)); end
else
  if pr, fprintf('\tError-order >= %i',min(4,p)+(~least_squares_appr)); end
end
if pr, fprintf('\n'); end

if (nargout>1) % Old syntax.
  w0 = w(1);
  w = w(2:end);
end


%%%%%%%%
function w=simple_calc(P)

n = size(P,2);

A = [ones(1,n); P; P.^2; P(1,:).*P(2,:)];
b = [0; 0;0; 2;2;0];
[U,S,V] = svd(A);
K = rank(S);
if (K<6) % rank K<6: there may be no solutions.
  dS = diag(S);
  w = V(:,1:K)*(diag(1./dS(1:K))*U(:,1:K)'*b);
  residual = A*w-b;
  if (norm(residual)>1e-6) % Underspecified configuration
    if (K>=4)
      if (K==5)
        % Remove x*y-eqn
        A_ = A(1:5,:);
        b_ = b(1:5);
      else % K==4
        % Remove x*y-eqn, and add the x^2 and y^2-eqns
        A_ = [A(1:3,:);sum(A(4:5,:),1)];
        b_ = [b(1:3);sum(b(4:5),1)];
      end
      [U,S,V] = svd(A_);
      dS = diag(S);
      w = V(:,1:K)*(diag(1./dS(1:K))*U(:,1:K)'*b_);
      residual2 = A*w-b;
      warning('IGMRF:underspecifiedConfiguration',...
              ['Rank deficient configuration.\n',...
               '  Residual_basic = %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n',...
               '  Residual_final = %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f'],...
              residual(1),residual(2),residual(3),residual(4),...
              residual(5),residual(6),...
              residual2(1),residual2(2),residual2(3),residual2(4),...
              residual2(5),residual2(6));
    else
      warning('IGMRF:underspecifiedConfiguration',...
              ['Rank deficient configuration.\n',...
               '  Residual = %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f'],...
              residual(1),residual(2),residual(3),residual(4), ...
              residual(5),residual(6));
    end
  end
else % rank K=6 (implies n>=6): the solution space has dimension n-K
  dS = diag(S);
  w = V(:,1:K)*(diag(1./dS)*U'*b);

  % If n==6, the solution is unique.
  % If n>=7 the solution space has dimension n-6.
  % w is here the solution with minimum w'*w,
  % which yields the minimum variance Laplace estimates, so don't proceed!
  if (n>=7) & (1==0)
    % Solution space: w + V(:,7:n)*t
    % Choose the weights with smallest difference:
    t = - ( (V(:,7:n)-ones(n,1)*mean(V(:,7:n))) \ ...
            (w-mean(w)) );
    w = w+V(:,7:n)*t;
  end
end




%%%%%%%%
function w=isotropic_error_calc(P)

n = size(P,2);

A = [ones(1,n); P; ...
     P(1,:).^2; P(1,:).*P(2,:); P(2,:).^2; ...
     P(1,:).^3; P(1,:).^2.*P(2,:); P(1,:).*P(2,:).^2; P(2,:).^3; ...
     P(1,:).^4; P(1,:).^3.*P(2,:); ...
     P(1,:).^2.*P(2,:).^2; ...
     P(1,:).*P(2,:).^3; P(2,:).^4];
A = [A, [zeros(1+2+3+4,1); 48;0;8;0;48]];
b = [0; 0;0; 2;0;2; 0;0;0;0; 0;0;0;0;0];
[U,S,V] = svd(A);
K = rank(S)
if (K<15) % rank K<15: there may be no solutions.
  dS = diag(S);
  w = V(:,1:K)*(diag(1./dS(1:K))*U(:,1:K)'*b);
  residual = A*w-b;
  if (norm(residual)>1e-6) % Underspecified configuration
    warning('IGMRF:underspecifiedConfiguration',...
            ['Rank deficient configuration.\n',...
             '         Residual = %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f'],...
            residual(1),residual(2),residual(3),residual(4), ...
            residual(5),residual(6));
  end
else % rank K=14 (implies d>=14): the solution space has dimension d-K
  dS = diag(S);
  w = V(:,1:14)*(diag(1./dS)*U'*b);
end

w0 = w(1);
w = w(2:end);
