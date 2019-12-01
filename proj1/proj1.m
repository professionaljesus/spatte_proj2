clear
close all
load('UStemp.mat')

%%%%%%%%%%%%%%%%%%OLS%%%%%%%%%%%%%%%%%%%
% Y(s) = X(s)*B + eta
% Residual noisy, eta ~sigma_e
% Do binned covariance estimation
% Re estimate the mean (B) with GLS

%%%%%%%%%%%%%%%%%%Kriging%%%%%%%%%%%%%%
% Y(s) = X(s)*B + eta
% ML is biased
% eta ~N(nu,sigma), [par,beta] = covest_ml, beta is regression term(?).

%Do OLS
%----need to check which variable are necessary.-----
X_reg = [ones(length(X),1) X];
beta = (X_reg'*X_reg)\(X_reg'*Y);
Y_pred = X_reg*beta;

%Get eta, should have zero mean and if no spatial dependence iid
%----need to evaluate spatial dependence-----
eta = Y - Y_pred;
sigma_eta = norm(eta)^2/(size(X,1)-size(X,2));
sigma_beta = (sigma_eta^2)*(X_reg'*X_reg)\eye(size(X_reg,2));

%Create distance matrix for known data stations
D = distance_matrix(X(:,1:2));

%Do binned covariance estimate, non-parametric
Dmax = max(max(D))*0.7;
Kmax = 45;
[rhat,s2hat,m,n,d] = covest_nonparametric(D,eta,Kmax,Dmax);

%Plot 1.covariance cloud with rhat in 2. binned covariance estimate 3. oberservation per bin.
plot(D,eta*eta','.k')
figure
plot(d,rhat,'-',0,s2hat,'o')
figure 
plot(m)
figure
%%
%----need to to bootstrap-----

%----Do WLS/LS estimation to obtian the parameter to the covariance function.

par_wls = covest_ls(rhat, s2hat, m, n, d, 'matern');

%Get covariance function for X_known
Sigma = matern_covariance(D,par_wls(1),par_wls(2),par_wls(3));

%Plot 1. covariance function 2. Above but now with the estimated covariance
%function
%----need to fix sigma plot like lecture Markov fields 1
plot(D,Sigma,'o')
figure 
plot(D,eta*eta','.k')
hold on 
plot(D,Sigma,'o')
figure
plot(d,rhat,'-',0,s2hat,'o')
hold on
plot(D,Sigma,'o')
figure

%GLS
%----needs to estimate relevant parameters & know what to do GLS on
Sigma = Sigma + s2hat^2*eye(length(X_reg));
beta_gls = (X_reg'*(Sigma\X_reg))\(X_reg'*(Sigma\Y));
Y_gls_reg = X_reg*beta_gls;
eta2 = Y-Y_gls_reg;
[rhat,s2hat,m,n,d] = covest_nonparametric(D,eta2,Kmax,Dmax);
par_gls = covest_ls(rhat, s2hat, m, n, d, 'matern');
Sigma = matern_covariance(D,par_gls(1),par_gls(2),par_gls(3));
Sigma = Sigma + s2hat^2*eye(length(X_reg));
beta_gls = (X_reg'*(Sigma\X_reg))\(X_reg'*(Sigma\Y));
Y_gls_reg = X_reg*beta_gls;
eta2 = Y-Y_gls_reg;

%Check with Validation
%1. Naive linear model
X_valdi_reg = [ones(length(X_valid),1) X_valid];
Y_pred = X_valdi_reg*beta;
wrong_naive = norm(Y_valid-Y_pred)^2;

%2. GLS
X_gls = [X(:,1:2); X_valid(:,1:2)];
D = distance_matrix(X_gls);
Sigma = matern_covariance(D,par_wls(1),par_wls(2),par_wls(3));
Sigma = Sigma + s2hat^2*eye(length(X_gls));

Sigma_kk = Sigma(1:500,1:500);
Sigma_uu = Sigma(501:600,501:600);
Sigma_ku = Sigma(1:500,501:600);
Sigma_uk = Sigma(501:600,1:500);


Y_pred = X_valdi_reg*beta_gls + Sigma(501:600,1:500)*(Sigma(1:500,1:500)\(Y-X_reg*beta_gls));
wrong_GLS = norm(Y_valid-Y_pred)^2;

var_gls = Sigma_uu - Sigma_uk*(Sigma_kk\Sigma_ku) + ... 
    (X_valdi_reg' - X_reg'*(Sigma_kk\Sigma_ku))'...
    *(((X_reg'*(Sigma_kk\X_reg)))...
    \(X_valdi_reg' - X_reg'*(Sigma_kk\Sigma_ku)));

%Reconstruct field
%(???) Use same Sigma_kk as above or it's the same with the new generated?
%(???) Need to add sigma_epsilon, prob not
%Make america great again
I_land = ~any(isnan(X_grid),2);
X_US_land = [X(:,1:2); X_grid(I_land,1:2)];
D = distance_matrix(X_US_land);
Sigma = matern_covariance(D,par_wls(1),par_wls(2),par_wls(3));
Sigma = Sigma + s2hat^2*eye(length(X_US_land));

X_grid_U = [ones(sum(I_land),1) X_grid(I_land,:)];
Y_grid = nan(sz_grid);
Y_grid(I_land) = X_grid_U*beta_gls+Sigma(501:end,1:500)*(Sigma(1:500,1:500)\(Y-X_reg*beta_gls));


imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid, ...
  'alphadata', reshape(I_land,sz_grid))
hold on
scatter(X(:,1), X(:,2), 20, Y, ...
  'filled','markeredgecolor','k')
scatter(X_valid(:,1), X_valid(:,2), 20, Y_valid, ...
  'filled','markeredgecolor','r')
axis xy tight; hold off; colorbar
title('GLS predictions')
figure
%Do Kriging
%---need to fix the uncertainty in the covariance estimation

D = distance_matrix(X(:,1:2));
[par_ml,beta_ml] = covest_ml(D,Y,'matern',[],X_reg);

%Kriging predictor, like with gsl
X_kriging = [X(:,1:2); X_valid(:,1:2)];
D = distance_matrix(X_kriging);
Sigma = matern_covariance(D,par_ml(1),par_ml(2),par_ml(3));
Sigma = Sigma + par_ml(4)^2*eye(length(X_kriging));

Sigma_kk = Sigma(1:500,1:500);
Sigma_uu = Sigma(501:600,501:600);
Sigma_ku = Sigma(1:500,501:600);
Sigma_uk = Sigma(501:600,1:500);

Y_pred = X_valdi_reg*beta_ml + Sigma_uk*(Sigma_kk\(Y-X_reg*beta_ml));
wrong_ML = norm(Y_valid-Y_pred)^2;

var_kriging = Sigma_uu - Sigma_uk*(Sigma_kk\Sigma_ku) + ... 
    (X_valdi_reg' - X_reg'*(Sigma_kk\Sigma_ku))'...
    *(((X_reg'*(Sigma_kk\X_reg)))...
    \(X_valdi_reg' - X_reg'*(Sigma_kk\Sigma_ku)));

%Field of USA_ml
D = distance_matrix(X_US_land);
Sigma = matern_covariance(D,par_ml(1),par_ml(2),par_ml(3));
Sigma = Sigma + par_ml(4)^2*eye(length(X_US_land));
Y_grid2 = nan(sz_grid);
Y_grid2(I_land) = X_grid_U*beta_ml+Sigma(501:end,1:500)*(Sigma(1:500,1:500)\(Y-X_reg*beta_ml));

imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid2, ...
  'alphadata', reshape(I_land,sz_grid))
hold on
scatter(X(:,1), X(:,2), 20, Y, ...
  'filled','markeredgecolor','k')
scatter(X_valid(:,1), X_valid(:,2), 20, Y_valid, ...
  'filled','markeredgecolor','r')
axis xy tight; hold off; colorbar
title('ML predictions')
figure

imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid2-Y_grid, ...
  'alphadata', reshape(I_land,sz_grid))
axis xy tight; colorbar
title('Diff')


