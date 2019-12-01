function y=fillholes(x)
% FILLHOLES Fill NaN-holes in an image.
%
%  y=fillholes(x)
%
% All NaN-values in x are replaced by the expectation
% of a 4-neighbour GMRF, conditional on the non-NaN-values.

% $Revision: 2952 $  $Date: 2006-09-22 18:45:42 +0200 (fre, 22 sep 2006) $
% Copyright (c) Finn Lindgren 2002-2003

[ny,nx,D] = size(x);

unknown = logical(any(isnan(x),3));
UK = unknown(:);

if ~any(UK)
  y = x;
  return
elseif all(UK)
  y = x;
  return
end

c = [0 1 0;1 0 1;0 1 0]/4;
edge = (conv2(unknown*1,c,'same')>0) & ~unknown;
E = edge(:);
tau2 = 1;
Q = gmrfprec(c,tau2,[ny,nx],(UK|E));

x = colstack(x);
y = x;
mu = mean(x(E,:),1);
muE = repmat(mu,[sum(E),1]);
muUK = repmat(mu,[sum(UK),1]);
y(UK,:) = muUK-Q(UK,UK)\(Q(UK,E)*(x(E,:)-muE));
y = icolstack(y,[ny,nx]);



function out=gmrfprec(c,tau2,sz,mask)
% GMRFPREC Compute precision matrices for GMRF:s
%
% Q=gmrfprec(c,tau2,sz,mask)
%
% c,tau2: GMRF-parameters, as defined by e.g. gmrfparamest.
% sz: The size of the lattice for which to compute the
%     covariance/precision matrix.
% mask: only compute the rows of Q for which mask>0

if (nargin<4), mask = []; end

np = sz(1)*sz(2); % Number of points.

out = sparse(np,np);
d1 = (size(c,1)-1)/2;
d2 = (size(c,1)-1)/2;
[I,J] = ndgrid(1:sz(1),1:sz(2));
I = I(:);
J = J(:);
if ~isempty(mask)
  I = I(mask(:));
  J = J(mask(:));
end
nij = length(I);
nc = sum(c(:)~=0);
II0 = zeros(nij,nc+1);
JJ0 = zeros(nij,nc+1);
II = zeros(nij,nc+1);
JJ = zeros(nij,nc+1);
value = zeros(nij,nc+1);
idx = 0;
cc = conv2(ones(sz(1),sz(2)),c,'same');
weight = sum(c(:))./(cc(:)+(cc(:)==0));
if ~isempty(mask)
  weight = weight(mask(:));
end
for i=-d1:d1
  for j=-d2:d2
    if (i==0) & (j==0)
      idx = idx+1;
      II0(:,idx) = I; JJ0(:,idx) = J;
      II(:,idx) = I;  JJ(:,idx) = J;
      ok(:,idx) = logical(1);
      value(:,idx) = 1/tau2;
    elseif (c(1+d1+i,1+d2+j)~=0)
      idx = idx+1;
      II0(:,idx) = I;  JJ0(:,idx) = J;
      II(:,idx) = I+i; JJ(:,idx) = J+j;
      ok(:,idx) = ((II(:,idx)>=1) & (II(:,idx)<=sz(1)) & ...
                   (JJ(:,idx)>=1) & (JJ(:,idx)<=sz(2)));
      value(:,idx) = -c(1+d1+i,1+d2+j)*weight/tau2;
    end
  end
end
out = sparse(II0(ok)+(JJ0(ok)-1)*sz(1),...
             II(ok)+(JJ(ok)-1)*sz(1),...
             value(ok),...
             np,np);
