function x=vec(X)
% VEC Vectorise a landmark matrix.
%
%  x=vec(X)
%
%  X: n landmark shapes, a p-by-m-by-n matrix
%  x: vectorised landmark shapes, a pm-by-n matrix

% $Id: vec.m 2952 2006-09-22 16:45:42Z johanl $

[p,m,n] = size(X);
x = reshape(X,[p*m,n]);
