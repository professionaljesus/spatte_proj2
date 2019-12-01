function x=reparametise(x0)
% REPARAMETISE Reparamatises a shape to obtain equidistant landmarks.
%
%  x=reparametise(x0)
%  x0 is a p-d-n matrix of landmarks and
%  size(x)=size(x0)

% $Id: reparameterise.m 4586 2012-10-08 16:18:33Z johanl $

x = zeros(size(x0));
n = size(x,1);
for j=1:size(x0,3)
  %calculate length between landmarks:
  L = diff(x0(:,:,j));
  L = sqrt(sum(L.^2,2));
  L = [L; sqrt( sum((x0(1,:,j)-x0(end,:,j)).^2) )];
  L = L/sum(L);
  Lc = cumsum(L);
  Lc = [-L(end);0;Lc;1+L(1)];

  xt = [x0(end,:,j);x0(:,:,j);x0(1:2,:,j)];
  
  ind = (0:n-1)/n;
  for i=1:size(x,2);
    x(:,i,j) = interp1(Lc,xt(:,i),ind,'spline');
  end;
end;
