function X=ivec(x,m)
% IVEC The inverse of vec, anti-vectorise landmarks.
%
%  X=ivec(x,m)
%
%  x: vectorised landmark shapes, a pm-by-n matrix
%  m: The landmark dimension.
%  X: n landmark shapes, a p-by-m-by-n matrix

% $Id: ivec.m 2952 2006-09-22 16:45:42Z johanl $

[pm,n] = size(x);
X = reshape(x,[pm/m,m,n]);
