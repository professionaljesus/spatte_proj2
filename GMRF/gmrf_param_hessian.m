function H=gmrf_param_hessian(gmrf_negloglike,theta0)
% GMRF_PARAM_HESSIAN  Calculates the hessian for a negated log-likelihood
%
% Example
%  H = gmrf_param_hessian(@(th) gmrf_negloglike(th,y,A,C,G,G2,B),...
%                         theta0)
%
% Number of function calls used:
%  1+2*dim+4*dim*(dim-1)/2 = 1+2*dim^2

% $Id: gmrf_param_hessian.m 4836 2014-12-10 11:09:32Z johanl $

delta = 1e-4;
dim = size(theta0,1);
H = zeros(dim,dim);
v0 = gmrf_negloglike(theta0);
for i=1:dim
  th = theta0;
  th(i) = theta0(i)-delta;
  v1 = gmrf_negloglike(th);
  th(i) = theta0(i)+delta;
  v2 = gmrf_negloglike(th);
  H(i,i) = (v1-2*v0+v2)/delta^2;
  
  for j=i+1:dim
    th = theta0;
    th(i) = theta0(i)-delta;
    th(j) = theta0(j)-delta;
    v11 = gmrf_negloglike(th);
    th(j) = theta0(j)+delta;
    v12 = gmrf_negloglike(th);
    th(i) = theta0(i)+delta;
    th(j) = theta0(j)-delta;
    v21 = gmrf_negloglike(th);
    th(j) = theta0(j)+delta;
    v22 = gmrf_negloglike(th);
    H(i,j) = (v22-v12-v21+v11)/(4*delta^2);
    H(j,i) = H(i,j);
  end
end
