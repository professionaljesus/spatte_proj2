function x=icolstack(y,sz);
% ICOLSTACK Invert column stacking of an image 
%
%  x=icolstack(y,sz)
%
%   sz: [m,n], size(x) = [sz,size(y,2)]
%       or
%       sz = size(x)
%
% SEE ALSO colstack

% $Id: icolstack.m 3318 2007-04-04 20:31:03Z finn $

x=reshape(y,[sz(1),sz(2),size(y,2)]);
