function Q = snake_neg_loglike(shape,shape_param,x,x_param)
% SNAKE_NEG_LOGLIKE Calculate negative log-likelihood for a snake.
%
% neg_logL = snake_neg_loglike(shape,shape_param)
% neg_logL = snake_neg_loglike(shape,shape_param,x,x_param)
%
% The first form computes "-log(likelihood for shape)"
% The second form also includes the data likelihood for the image x.
%
% shape: p-by-2-matrix with landmarks controlling a spline curve.
% shape_param: A struct with fields alpha, beta, and subdiv;
%   shape_param.alpha  : Squared derivative coefficient, default 0.005.
%   shape_param.beta   : Squared 2nd derivative coefficient, default 0.00001.
%   shape_param.subdiv : Number of spline subdivision levels
%                        (see simplespline, default=3)
% 
% x: A grayscale image.
% x_param: A struct with fields mu and sigma2
%   x_param.mu     : mu(1) is the image expectation outside the shape.
%                    mu(2) is the expectation inside the shape.
%   x_param.sigma2 : sigma2(1) is the variance outside the shape.
%                    sigma2(2) is the variance inside the shape.
%   Defaults: mu = [0,1], sigma2 = [1,1];

% Copyright (c) Finn Lindgren 2001-2003,2004
% $Id: snake_neg_loglike.m 2957 2006-09-25 07:20:23Z johanl $

if (nargin<2), shape_param = []; end
if (nargin<4), x_param = []; end

shape_param = defstruct(shape_param,struct(...
    'alpha',1e-3,...
    'beta',1e-3,...
    'subdiv',3));

shape_spline = simplespline(shape,shape_param.subdiv);
Q = Q_omega(shape_spline,shape_param);

if (nargin>=3)
  x_param = defstruct(x_param,struct(...
      'mu',[0,1],...
      'sigma2',[1,1],...
      'sz',size(x)));
  [mu,sigma2] = expect(shape_spline,x_param);
  Q = Q + Q_x(x,mu,sigma2);
end



function Q = Q_omega(shape_spline,shape_param)

Q = shape_param.alpha * size(shape_spline,1) * ...
    sum(sum((shape_spline - shape_spline([end,1:end-1],:)).^2,2),1);

Q = Q + shape_param.beta * size(shape_spline,1)^3 * ...
    sum(sum((shape_spline([2:end,1],:) - ...
             2*shape_spline + ...
             shape_spline([end,1:end-1],:)).^2,2),1);

Q = Q/2;



function Q = Q_x(x,mu,sigma2)

Q = sum(sum(...
    (x-mu).^2./(2*sigma2) +...
    0.5*log(2*pi) + 0.5*log(sigma2) ));


function [mu,sigma2]=expect(shape_spline,x_param);

I = indicshape(shape_spline,x_param.sz);
mu = (1-I)*x_param.mu(1) + I*x_param.mu(2);
sigma2 = (1-I)*x_param.sigma2(1) + I*x_param.sigma2(2);
