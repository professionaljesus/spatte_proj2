function Z0=procrustes_mean(Z,method,reflections)
% PROCRUSTES_MEAN Estimate the Procrustes mean.
%
%  Z0=procrustes_mean(Z,method,reflections)
%
%  Z should be preshapes, either Helmertised or centred.
%  method: 'full', full Procrustes mean. (default)
%          'partial', partial Procrustes mean.
%  reflections: 1, reflection is allowed. (default)
%               0, only pure rotations are allowed.

% Copyright (c) Finn Lindgren, 2005
% $Id: procrustes_mean.m 3233 2007-02-07 23:20:43Z finn $

if (nargin<2), method = []; end
if (nargin<3), reflections = []; end
if (isempty(method)), method = 'full'; end
if (isempty(reflections)), reflections = 1; end

[p,m,n] = size(Z);

% Use the first shape as initial estimate:
Z0 = Z(:,:,1);
beta_mean_0 = -1;
beta_mean_1 = 0;
while (abs(beta_mean_1-beta_mean_0)>1e-10)
  % Align the shapes:
  [ZP,beta] = procrustes_align(Z,Z0,method,reflections);
  % Prepare the stopping criterion:
  beta_mean_0 = beta_mean_1;
  beta_mean_1 = mean(beta);
  % Compute the new mean:
  Z0 = mean(ZP,3);
  % Normalise the size:
  Z0 = Z0/sqrt(sum(Z0(:).^2));
end
