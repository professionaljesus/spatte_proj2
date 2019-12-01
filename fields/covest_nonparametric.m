function [rhat,s2hat,m,n,d,varioest] = covest_nonparametric(D,z,Kmax,Dmax)
% covest_nonparametric  nonparametric covariance estimator
%
% Alt 1: [rhat,s2hat,m,n,d]=covest_nonparametric(D,z,Kmax,Dmax)
% Alt 2: [rhat,s2hat,m,n,d]=covest_nonparametric(U,z,Kmax,Dmax)
%
% Optional:
%   [...,varioest]=covest_nonparametric(...)
% also calculates the semi-variogram
%
%
% Alt 1 is slow, and requires O(n^2) memory
% Alt 2 is quicker, and requires only O(n) memory
%
% rhat,s2hat = estimated covariance and data variance
%              bins with no observations will be = NaN
% m = the number of observation pairs used for each bin
% n = the number of observations used for the s2hat estimation
% d = the distances for rhat
%
% D = the distance matrix for the data (n-by-n, for Alt 1)
% U = the coordinate matrix (n-by-dimension, for Alt 2)
% z = the residuals, z = y - mu, for fixed mu
% Kmax = the number of bins is Kmax+1
%
% Dmax = the maximum distance; the last bin ends at Dmax
%        Default = max(D(:))+0.001 (to ensure that all observations are
%        used)
%        Dmax is optional for Alt 1, but required for Alt 2.
%
% Example plot:
%   plot(d,rhat,'-',0,s2hat,'o')
% And for the variogram
%   plot(d,varioest,'-')

% $Id: covest_nonparametric.m 5090 2017-11-05 19:35:05Z johanl $

if (nargin<4), Dmax = []; end

if (size(D,1)~=size(z,1))
  error(['The number of coordinates, %d%, does not match the ' ...
		'number of data points, %d%.'],size(D,1),size(z,1))
end
if (size(D,1)==size(D,2))
  altmethod = 1;
  if isempty(Dmax), Dmax = max(D(:))+0.001; end
else
  altmethod = 2;
  if isempty(Dmax), error('Dmax must be supplied'), end
end

h = Dmax/(Kmax+1);
% bin-boundaries
d = h*(0:Kmax+1);
m = zeros(1,Kmax+1);
rhat = zeros(1,Kmax+1);
if (nargout>=6)
  varioest = zeros(1,Kmax+1);
end
n = length(z);
s2hat = z'*z/n;

if (altmethod==1)
  for k=0:Kmax
    ind = (d(k+1)<=D) & (D<d(k+2));
    [I,J] = find(ind);
    ok = (I<J); % Remove i==j, and the duplicate pairs (i,j)<-->(j,i)
    m(k+1) = sum(ok);
    rhat(k+1) = z(I(ok))'*z(J(ok))/m(k+1);
    if (nargout>=6)
      varioest(k+1) = 0.5*sum((z(I(ok))-z(J(ok))).^2)/m(k+1);
    end
  end
else % (altmethod==2)
  for i=1:size(D,1)
    D_ = sqrt(sum((D(i+1:end,:)-ones(n-i,1)*D(i,:)).^2,2));
    k = floor(D_/Dmax*(Kmax+1));
    ok = find(k<=Kmax);
    nok = length(ok);
    m = m+sparse(ones(nok,1),k(ok)+1,ones(nok,1),1,Kmax+1);
    rhat = rhat+sparse(ones(nok,1),k(ok)+1,z(i)*z(i+ok),...
                       1,Kmax+1);
    if (nargout>=6)
      varioest = varioest +...
          sparse(ones(nok,1),k(ok)+1,(z(i)-z(i+ok)).^2,...
                 1,Kmax+1);
    end
  end
  rhat = rhat./m;
  if (nargout>=6)
    varioest = 0.5*(varioest./m);
  end
end

d = d(1:end-1);
