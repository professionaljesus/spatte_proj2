function [ZP,beta,Gamma]=procrustes_align(Z,Z0,method,reflections)
% PROCRUSTES_ALIGN Perform Procrustes alignment of landmark data.
%   Aligns landmark data to a given template shape.
%
%  [ZP,beta,Gamma]=procrustes_align(Z,Z0)
%  [ZP,beta,Gamma]=procrustes_align(Z,Z0,method)
%  [ZP,beta,Gamma]=procrustes_align(Z,Z0,method,reflections)
%
%  Z:  n preshapes, either Helmertised or centred.
%  Z0: the preshape to align the Z-shapes to.
%  method: 'full', computes full Procrustes alignments. (default)
%          'partial', computes partial Procrustes alignments. 
%  reflections: 1, reflection-space is used. (default)
%               0, only rotations are allowed.
%
%  ZP: n aligned preshapes
%  beta: The scaling factors for full Procrustes alignment, 1-by-n
%  Gamma: The rotation matrices, cell array, 1-by-n

% Copyright (c) 2002-2004 by Finn Lindgren
% $Id: procrustes_align.m 3725 2008-12-02 08:57:38Z johanl $

if (nargin<3), method = []; end
if (nargin<4), reflections = []; end
if isempty(method), method = 'full'; end
if isempty(reflections), reflections = 1; end

[p,m,n] = size(Z); % p = number of landmarks or one less.

ZP = Z;
for k=1:n
  [U,S,V] = svd(Z(:,:,k)'*Z0);
  Gamma{k} = U*V';
  if (reflections)
    beta(1,k) = sum(diag(S));
  else
    if (det(Gamma{k})<0)
      Gamma{k} = [U(:,1:end-1),-U(:,end)]*V';
      beta(1,k) = sum(diag(S))-2*S(end,end);
    else
      beta(1,k) = sum(diag(S));
    end
  end
  switch (method)
   case 'full',
    ZP(:,:,k) = Z(:,:,k)*Gamma{k}*beta(k);
   case 'partial'
    ZP(:,:,k) = Z(:,:,k)*Gamma{k};
  end
end
