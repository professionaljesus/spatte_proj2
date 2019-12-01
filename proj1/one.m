% First use matern_covariance to create a Sigma-covariance matrix.
% and set mu=(a constant mean).
[u1,u2] = ndgrid(1:50,1:60);
D = distance_matrix([u1(:), u2(:)]);
Sigma = matern_covariance(D,1,1,1);
mu = 0.1;
N = 3000;
sz = [50 60];
R = chol(Sigma); % Calculate the Cholesky factorisation
eta = mu+R'*randn(N,1); % Simulate a sample
eta_image = reshape(eta,sz); %reshape the column to an image

imagesc(eta_image)