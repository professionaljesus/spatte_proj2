function y=colstack(x);
% COLSTACK Stack columns of an image 
%
%  y=colstack(x)
%
% SEE ALSO icolstack

% $Id: colstack.m 3318 2007-04-04 20:31:03Z finn $

y=reshape(x,[size(x,1)*size(x,2),size(x,3)]);
