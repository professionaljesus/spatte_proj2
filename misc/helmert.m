function H=helmert(k,sub)
% HELMERT Construct a Helmert (sub-)matrix.
% 
%  H = helmert(p)  or  H = helmert(p,'sub')
%   gives a Helmert sub-matrix with  p-1  rows and  p  columns.
%  H = helmert(p,'full')
%   gives a full Helmert matrix with  p  rows and columns

% $Id: helmert.m 2952 2006-09-22 16:45:42Z johanl $

if (nargin<2), sub = []; end
if isempty(sub), sub = 'sub'; end

switch sub
 case 'sub',
  H = -gallery('orthog',k,4);
  H(1,:) = [];
 case 'full',
  H = -gallery('orthog',k,4);
  H(1,:) = H(1,:)*sign(H(1,1));
 otherwise,
  error('Second parameter must be [], ''sub'', or ''full''.')
end
