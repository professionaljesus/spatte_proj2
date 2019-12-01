function Z=preshape(X,centred)
% PRESHAPE Compute the preshapes of landmark objects
% 
%  Z = preshape(X)   gives Helmertised preshapes,
%                    Use H'*Z(:,:,k) to obtain icons,
%                    where H is computed by the helmert function.
%  Z = preshape(X,1) gives centred preshapes
%                    Z(:,:,k) can be used as icons.
%
%  X should be a p-by-m-by-n matrix, where each X(:,:,k) is a shape
%  in m dimensions, with p landmarks.

% $Id: preshape.m 2957 2006-09-25 07:20:23Z johanl $

if (nargin<2), centred = []; end
if isempty(centred), centred = 0; end

[p,m,n] = size(X);

A = helmert(p);
if (centred)
  A = A'*A;
  Z = zeros(p,m,n);
else
  Z = zeros(p-1,m,n);
end
for k=1:size(X,3)
  tmp = A*X(:,:,k);
  Z(:,:,k) = tmp/sqrt(sum(tmp(:).^2));
end
