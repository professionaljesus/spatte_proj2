function distance=procrustes_dist(beta)
% PROCRUSTES_DIST Compute Procrustes distances.
%
%  distance=procrustes_dist(beta)
%
%  beta should be the scaling factors computed in the
%  full Procrustes alignment, as output from  procrustes_align
%  distance(1,:) are the full Procrustes distances
%  distance(2,:) are the partial Procrustes distances
%  distance(3,:) are the Procrustes distances

% $Id: procrustes_dist.m 2957 2006-09-25 07:20:23Z johanl $

n = length(beta);

distance = zeros(3,n);
distance(1,:) = sqrt(1-beta.^2);
distance(2,:) = sqrt(2-2*beta);
distance(3,:) = acos(beta);
