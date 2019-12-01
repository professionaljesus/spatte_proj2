function Q=gmrfprec(sz,q)
% GMRFPREC Constructs a precision matrix for a GMRF on a regular grid
%
%   Q = gmrfprec(sz,q)
%
%   sz : The grid/image size,  sz = [m,n] for an m-by-n image
%   q  : The precision specification.  Example:
%        q = [ 0  -1   0;...
%             -1   5  -1;...
%              0  -1   0];
%   Q  : A sparse precision matrix
%
%   No special care is taken on the grid borders.
%   The resulting model is equivalent to saying that the field is equal
%   to its expectation outside the grid.
%
% See also: igmrfprec

% $Id: gmrfprec.m 4836 2014-12-10 11:09:32Z johanl $

n = prod(sz);
qsz = size(q);
qn = (qsz-1)/2;

II_I = [];
II_J = [];
JJ_I = [];
JJ_J = [];
KK = [];
[I,J] = ndgrid(1:sz(1),1:sz(2));
I = I(:); J = J(:);
for a=-qn(1):qn(1)
  for b=-qn(2):qn(2)
    if (q(qn(1)+1+a,qn(2)+1+b) ~= 0)
      II_I = [II_I;I];
      II_J = [II_J;J];
      JJ_I = [JJ_I;I+a];
      JJ_J = [JJ_J;J+b];
      KK = [KK;...
            ones(length(I),1)*q(qn(1)+1+a,qn(2)+1+b)];
    end
  end
end
II = II_I+sz(1)*(II_J-1);
JJ = JJ_I+sz(1)*(JJ_J-1);
ok = (II_I>=1) & (II_I<=sz(1)) & (II_J>=1) & (II_J<=sz(2)) & ...
     (JJ_I>=1) & (JJ_I<=sz(1)) & (JJ_J>=1) & (JJ_J<=sz(2));
II(~ok) = [];
JJ(~ok) = [];
KK(~ok) = [];
Q = sparse(II,JJ,KK,n,n);
