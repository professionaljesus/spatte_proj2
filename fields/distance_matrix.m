function D=distance_matrix(u,v)
% distance_matrix  calculates the distance matrix for one or two sets of locations
%
%   D = distance_matrix(u)
%   D = distance_matrix(u, v)
%
% D(i,j) = sqrt(sum((u(j,:)-u(i,:)).^2)) = || u_j - u_i ||
% u = N-by-d matrix, where N is the number of locations,
%     and d is the space dimension
%
% To compute the cross-distances between to sets of locations use:
%   D = distance_matrix(u, v)
%
% D(i,j) = sqrt(sum((u(j,:)-v(i,:)).^2)) = || u_j - v_i ||
% u = N1-by-d matrix, where N is the number of locations,
%     and d is the space dimension
% v = N2-by-d matrix, where N is the number of locations,
%     and d is the space dimension
%
% Example:, for a 50-by-60 image
%   [u1,u2] = ndgrid(1:50,1:60);
%   D = distance_matrix([u1(:), u2(:)]);
%
%   u = [u1(:), u2(:)];
%   crossD = distance_matrix(u(1:500,:), u(501:1000,:));
%   max(max( abs(crossD - D(1:500,501:1000)) ))


% $Id: distance_matrix.m 4984 2016-11-06 20:01:13Z johanl $

if nargin<2 || isempty(v), v=u; end

if size(u,2)~=size(v,2)
  error('second dimension of u and v must match')
end

D = zeros(size(u,1), size(v,1));
for k=1:size(u,2)
  D = D + bsxfun(@minus,u(:,k),v(:,k)').^2;
end
D = sqrt(D);
