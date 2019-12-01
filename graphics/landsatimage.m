function y=landsatimage(x,bands,outliers)
% LANDSATIMAGE Make an RGB image matrix from LANDSAT data.
%              Outliers are removed, and the intensities are rescaled.
%
%  y=landsatimage(x)
%  y=landsatimage(x,bands)
%
%  x: the LANDSAT image, m-n-7 matrix. Rescaled or unscaled.
%  bands: The spectral bands to use, 3-element vector.
%         Default = [3,2,1] (R,G,B)
%                   [4,3,2] highlights trees (infrared) in red,
%  y: Transformed image values, m-n-3 matrix, suitable for "image".
%
% Example:
%  image(landsatimage(lanread('rio.lan')));

% Requires colstack.

% Copyright (c) 2002 Finn Lindgren
% $Revision: 2952 $  $Date: 2006-09-22 18:45:42 +0200 (fre, 22 sep 2006) $

if nargin<2
  bands = [];
end
if nargin<3
  outliers = [];
end
if isempty(bands)
  bands = [3,2,1];
end
if isempty(outliers)
  outliers = 0.03;
end

if ~isa(x,'double')
  x = double(x);
end
x = x(:,:,bands);

[m,n,d] = size(x);
sortx = sort(colstack(x),1);
low_quantiles = sortx(max(1,floor(m*n*outliers)),:);
high_quantiles = sortx(min(m*n,ceil(m*n*(1-outliers))),:);
low_quantiles = repmat(reshape(low_quantiles,[1,1,d]),[m,n,1]);
high_quantiles = repmat(reshape(high_quantiles,[1,1,d]),[m,n,1]);
y = max(0,min(1,(x-low_quantiles)./(high_quantiles-low_quantiles)));
