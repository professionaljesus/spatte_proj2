function B = sphereHarmonics(x,y,z,order)
% SPHEREHARMONICS Creates spherical harmonic functions.
%
%  B = sphereHarmonics(x,y,z,order)
%  B = sphereHarmonics(P,order)
%  or
%  B = sphereHarmonics(long,lat,order)
%
%  Input should be column vector or matrices.
%
%  Constructs harmonics up to order for the input points.
%  P = [x,y,z] coordinates are assumed to be on a circle centred at the 
%  origin. 
%
%  long,lat coordinates are assumed to be in degrees.

% $Id: sphereHarmonics.m 4586 2012-10-08 16:18:33Z johanl $

if nargin==3 %long,lat
  cos_theta = cos(pi/2-y/180*pi);
  phi = x/180*pi;
  order = z;
else %points
  if nargin==4
    P = [x y z];
    clear x y z;
  else
    P = x;
    order = y;
  end
  %calculate radie
  R = sqrt(sum(P(1,:).^2,2));
  cos_theta = P(:,3)/R;
  phi = atan2(P(:,2),P(:,1));
end

N = length(phi);
B = zeros(N,(order+1)^2);

%first construct legandre polynomials by recursion
B(:,1) = 1;
if order==0
  return;
end
B(:,3) = cos_theta;
for k=2:order
  B(:,k.^2+k+1) = ( (2*k-1)*B(:,3).*B(:,k.^2-k+1) - (k-1)*B(:,k.^2-3*k+3) )/k;
end
%second construct associated legendre functions with m=l
B(:,4) = sqrt((1-cos_theta.^2));
for k=2:order
  B(:,k.^2+2*k+1) = prod((2*k-1):-2:1)*B(:,4).^k;
end
%third construct remaining functions for 0<m<l using recursion.
for k=2:order
  for l=1:k-1
    B(:,k.^2+k+1+l) = ((k+l-1)*B(:,k.^2-k+1+(l-1)) - (k-(l-1))*B(:,3).*B(:,k.^2+k+1+(l-1)))./B(:,4);
  end
end
%fourth multiply with sin or cos and scale
for k=1:order
  B(:,k.^2+k+1) = sqrt(2*k+1)*B(:,k.^2+k+1);
  for l=1:k
    MULT = sqrt((2*k+1)*factorial(k-l)/factorial(k+l));
    B(:,k.^2+k+1-l) = MULT * B(:,k.^2+k+1+l).*sin(l*phi);
    B(:,k.^2+k+1+l) = MULT * B(:,k.^2+k+1+l).*cos(l*phi);
  end
end
