function theta = normmix_initial(x, K, gt, varDiag, D)
% NORMMIX_INITIAL Computes initial estimates for K classes
%
% theta = normmix_initial(x, K, gt)
% theta = normmix_initial(x, K, gt, varDiag, D)
%
% x: n-by-d matrix
% K: number of classes
% gt: n-by-1 matrix with ground truth data (1..K indicating apriori known
%     class for a pixel, 0 indicating unknown class).
% varDiag: diagonal elements to add to each estimated covariance matrix.
%          Defaults to varDiag = var(x)/K.
% D: Threshold for considering pixels as belonging to a class with ground
%    truth data.
%    Defaults to D = chi2inv(.5, size(x,2))
%
% The function first computes estimates of mean and covariance matrices for
% pixels that have ground-truth. If K is larger than number of classes with
% specified ground truth the function then:
% 1) excludes pixels with values similar to these proto-classes, i.e. pixels
%    for which:  
%      (x-mu_k)' * inv(Sigma_k) * (x-mu_k) < D
% 2) Runs Kmeans to group the remaining pixels and estimates mean and
%    covariance matrices for each of the groups
%
% Returns theta with
% theta{k}.mu: 1-by-d matrix, class expected value.
% theta{k}.Sigma: d-by-d matrix, class covariance.

% $Id: normmix_initial.m 4792 2014-11-27 14:30:32Z johanl $

if nargin<4, varDiag=[]; end
if nargin<5, D=[]; end

if isempty(varDiag), varDiag = var(x)/K; end
if isempty(D), D = chi2inv(.5, size(x,2)); end

if isscalar(varDiag)
    varDiag = eye(size(x,2))*varDiag;
elseif length(varDiag)==size(x,2)
    varDiag = diag(varDiag);
else
    error('Bad size for varDiag')
end
theta = cell(1,K);

%indicator for pixels falling within a gt-specified class.
I = false(size(x,1),1);

%number of classes not specified through gt
Krest = K;

%compute estimates for each class based on gt.
for k=1:K
  if sum(gt==k)~=0
    theta{k}.mu = mean( x(gt==k,:) );
    theta{k}.Sigma = cov( x(gt==k,:) ) + varDiag;
    %reduce unknown classes by one
    Krest = Krest-1;
    %find pixels in this class.
    y = x-repmat(theta{k}.mu,[size(x,1) 1]);
    y = sum((y/theta{k}.Sigma).*y,2);
    Ik = y<D;
    I = I | Ik;
  end
end
if Krest==0
    %all classes estimated, DONE!
    return
end
%classify remaining pixels
y = x(~I,:);
cl = kmeans(y, Krest);
%estimate initial parameters for remaining classes based on Kmeans estimate
kInd = 1;
for k=1:K
  if sum(gt==k)==0
    theta{k}.mu = mean( y(cl==kInd,:) );
    theta{k}.Sigma = cov( y(cl==kInd,:) ) + varDiag;
    kInd = kInd+1;
  end
end
