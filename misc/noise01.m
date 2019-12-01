function y=noise01(x,p_error)
% NOISE01 Modify binary image with noise.
%
%  y=noise01(x,p_error);
%
%  x: binary image, x_u is 0 or 1.
%  y: binary image, y_u is 0 or 1.
%
%  P(y_u = x_u)  = 1-p_error
%  P(y_u ~= x_u) = p_error

% $Revision: 2952 $

y=(((x*2-1).*(1-2*(rand(size(x))<p_error)))>0);

