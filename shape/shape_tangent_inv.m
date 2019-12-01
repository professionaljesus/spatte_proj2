function Z=shape_tangent_inv(v,Z0,method)
% SHAPE_TANGENT_INV Compute pre-shapes from shape tangent coordinates.
%
%  Z=shape_tangent_inv(v,Z0,method)
%
%  v should be vectorised tangent coordinates
%  Z0 should be the tangent plane pole preshape
%  Z will be the back-projected preshapes, either to the
%  preshape sphere (method='partial', default) or the
%  full Procrustes sphere (method='full')

% $Id: shape_tangent_inv.m 2957 2006-09-25 07:20:23Z johanl $

if (nargin<3), method = []; end
if (isempty(method)), method = 'partial'; end

[p,m] = size(Z0); % p = number of landmarks or one less.
n = size(v,2);

vecZ0 = vec(Z0);
if strcmp(method,'partial')
  z = vecZ0*sqrt(1-sum(v.^2,1))+v;
elseif strcmp(method,'full')
  z = vecZ0*(0.5+sqrt(0.25-sum(v.^2,1)))+v;
end
Z = ivec(z,m);
