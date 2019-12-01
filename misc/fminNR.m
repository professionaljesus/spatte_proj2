function [x0, fval, exitflag] = fminNR(fun, x0, options)
% FMINNR Function minimisation using Newton-Raphson
%
%  [x,fval,exitflag] = fminNR(fun, x0, options)
%
%  fun: Function to minimize, should be of the form 
%       [f, df, d2f] = fun(x)
%       i.e. returning function value, derivative, and Hessian for each value x.
%  x0: initial values for the optimisation.

% $Id: fminNR.m 4955 2016-04-12 18:32:17Z johanl $

%options
tol_fun = 1e-7;
tol_X = 1e-7;
tol_df = 1e-7;
maxIter = 1e2;

%compute funtion at x0
[fval, df, d2f] = fun(x0);

%convergence diagnostics
done = false;
exitflag = 0;
i = 1;

while ~done && exitflag==0
  %save old values from last step
  x_old = x0;
  f_old = fval;
  df_old = df;
  
  %compute search direction
  [R, d2f_OK] = chol(d2f); 
  if d2f_OK==0
    p = -R\(df(:)'/R)';
  else
    keyboard
    p = df(:);
  end
  %compute new x-value
  x0 = x_old + reshape(p,size(x0));
  
  %compute new values (assuming standard NR step)
  [fval, df, d2f] = fun(x0);
  
  %compute p'*df for old and new df
  p_df_new = p'*df(:);
  p_df_old = p'*df_old(:);

  if ~checkWolfe(fval, f_old, p_df_new, p_df_old, 1) && norm(df)>tol_df
    %if not do a line search (maybe needed repeatedly)
    k=0;
    while (k==0 || fval>f_old) && k<10
      %shrink upper bound as 2^-k
      alpha = fminbnd(@(t) fun(x_old + t*reshape(p,size(x0))),0,2^(-k));
      %try new point
      x0 = x_old + alpha*reshape(p,size(x0));
      [fval, df, d2f] = fun(x0);
      k=k+1;
    end
    p_df_new = p'*df(:);
    if ~checkWolfe(fval, f_old, p_df_new, p_df_old, 1) && norm(df)>tol_df
      %still not converged -> flag as error
      exitflag = -1;
    end
  end

  if exitflag==0
    %have computations converged
    done = (abs(fval-f_old)/(abs(fval)+tol_fun) < tol_fun) || ...
          (norm(x0-x_old)/(norm(x0)+tol_X) < tol_X) ||  ...
          (norm(df) < tol_df);
     i = i+1;
    exitflag = i>maxIter;
  end
end

function OK = checkWolfe(f_new, f_old, p_df_new, p_df_old, alpha)
%check if the Wolfe conditions are fulfilled
c1 = 1e-4;
c2 = 0.9; %Newton-method
%1) f_new < f_old + c1*alpha * p'*df_old
OK = f_new < f_old + c1*alpha*p_df_old;
%2) 
OK = OK && abs(p_df_new) < c2*abs(p_df_old);