function v=shape_tangent(Z,Z0,method,reflections)
% SHAPE_TANGENT Compute vectorised shape space tangent coordinates.
%   Aligns the preshapes in Z to Z0, and calculates tangent coordinates.
%
%  v=shape_tangent(Z,Z0,method)
%
%  Z should be preshapes, either Helmertised or centred.
%  Z0 should be the tangent plane pole preshape
%  Z will be Procrustes aligned to Z0 before projection onto
%  the tangent plane at Z0, generating the vectorised tangent
%  coordinates v.
%  method: 'partial', partial Procrustes tangent coordinates. (default)
%          'full', full Procrustes tangents coordinates.

% $Id: shape_tangent.m 3724 2008-12-02 08:55:57Z johanl $

if (nargin<3), method = []; end
if (nargin<4), reflections = []; end
if (isempty(method)), method = 'partial'; end
if (isempty(reflections)), reflections = 1; end

[p,m,n] = size(Z); % p = the number of landmarks or one less.

ZP = procrustes_align(Z,Z0,method,reflections);
vecZ0 = vec(Z0);
v = (eye(length(vecZ0))-vecZ0*vecZ0')*vec(ZP);
