function marks2=simplespline(marks0,levels)
% SIMPLESPLINE Compute a spline interpolation of a sequence of landmarks.
%   Method: Subdivide a closed curve using the 4-point method.
%
%  marks1=simplespline(marks0)
%  marks1=simplespline(marks0,levels)
%
%  marks0,marks1: (p-by-2)-matrices of landmark coordinates
%  levels: the number of subdivision levels, default=3.

% Copyright (c) Finn Lindgren 2001, 2004
% $Id: simplespline.m 2957 2006-09-25 07:20:23Z johanl $

if (nargin<2), levels = []; end
if (isempty(levels)), levels = 3; end

[p,d] = size(marks0);
marks2 = marks0;
for level=1:levels
  marks1 = marks2;
  marks2 = zeros(2*p,d);
  marks2(1:2:2*p-1,:) = marks1;
  marks2(2:2:2*p,:) = conv2(marks1([end,1:end,1:2],:),...
			 [-1;9;9;-1]/16,'valid');
  p = 2*p;
end
