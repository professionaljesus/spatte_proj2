% Load data
load fmri.mat

Xa = [X(:,1:2) sum(X(:,3:end),2)];

Y = img(:);
sz = size(img);
A = kron(Xa, speye(sz(1)*sz(2)));
%Gibb loop


%Posterior expectation/variance
% mean(beta) ?? std(beta)?????


%Significant pixel
