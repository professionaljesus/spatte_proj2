clear
close all
load('UStemp.mat')

%%%%%%%%%%%%%%%%%%OLS%%%%%%%%%%%%%%%%%%%
% Y(s) = X(s)*B + eta
% Residual noisy, eta ~sigma_e
% Do binned covariance estimation
% Re estimate the mean (B) with GLS

%%%%%%%%%%%%%%%%%%Kriging%%%%%%%%%%%%%%
% Y(s) = X(s)*B + eta + epsilon
% ML is biased
% eta ~N(nu,sigma), [par,beta] = covest_ml, beta is regression term(?).

%Do OLS
%----need to check which variable are necessary.-----
X_reg = [ones(length(X),1) X X(:,1).^2 X(:,2).^2 X(:,4).^2 X(:,5).^2];
X_valid_reg = [ones(length(X_valid),1) X_valid X_valid(:,1).^2 X_valid(:,2).^2 X_valid(:,4).^2 X_valid(:,5).^2];
beta_ols = (X_reg'*X_reg)\(X_reg'*Y);
Y_pred = X_reg*beta_ols;

%Get eta, should have zero mean and if no spatial dependence iid
%----need to evaluate spatial dependence-----
eta = Y - Y_pred;
sigma_eta = norm(eta)^2/(size(X_reg,1)-size(X_reg,2));
I_land = ~any(isnan(X_grid),2);

temp = X_grid(I_land,:);
sigma_beta = (sigma_eta)*(temp'*temp)\eye(size(temp,2));

Vmu = sum((temp*sigma_beta).*temp,2);
var_ols = sigma_eta +  Vmu;

Y_valid_pred = X_valid_reg*beta_ols;

normOLS = (norm(Y_valid_pred - Y_valid)^2)/100;

I_land = ~any(isnan(X_grid),2);
X_US_land = [X(:,1:2); X_grid(I_land,1:2)];

X_grid_U = [ones(sum(I_land),1) X_grid(I_land,:) X_grid(I_land,1).^2 X_grid(I_land,2).^2 X_grid(I_land,4).^2 X_grid(I_land,5).^2];
Y_grid = nan(sz_grid);
Y_grid_low = Y_grid;
Y_grid_hi = Y_grid;
Y_grid(I_land) = X_grid_U*beta_ols;
Y_grid_low(I_land) = X_grid_U*beta_ols - 1.96*sqrt(var_ols);
Y_grid_hi(I_land) = X_grid_U*beta_ols + 1.96*sqrt(var_ols);


figure
subplot(221)
imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid, ...
  'alphadata', reshape(I_land,sz_grid))
hold on
scatter(X(:,1), X(:,2), 20, Y, ...
  'filled','markeredgecolor','k')
scatter(X_valid(:,1), X_valid(:,2), 20, Y_valid, ...
  'filled','markeredgecolor','r')
caxis([-30 30]);
hold off; axis xy tight; colorbar
title('OLS predictions')

subplot(222)
imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid, ...
  'alphadata', reshape(I_land,sz_grid))
caxis([-30 30]);
hold off; axis xy tight; colorbar
title('OLS reconstruction')

subplot(223)
imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid_low, ...
  'alphadata', reshape(I_land,sz_grid))
caxis([-30 30]);

axis xy tight; colorbar
title('OLS low')

subplot(224)
imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid_hi, ...
  'alphadata', reshape(I_land,sz_grid))
caxis([-30 30]);

axis xy tight; colorbar
title('OLS hi')
%%
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

%Do Kriging
%---need to fix the uncertainty in the covariance estimation

D = distance_matrix(X(:,1:2));
[par_ml,beta_ml] = covest_ml(D,Y,'matern',[],X_reg);

%Kriging predictor, like with gsl
X_kriging = [X(:,1:2); X_valid(:,1:2)];
D = distance_matrix(X_kriging);
Sigma = matern_covariance(D,par_ml(1),par_ml(2),par_ml(3));
Sigma = Sigma + par_ml(4)*eye(length(X_kriging));
Sigma_kk = Sigma(1:500,1:500);
Sigma_uu = Sigma(501:600,501:600);
Sigma_ku = Sigma(1:500,501:600);
Sigma_uk = Sigma(501:600,1:500);

Y_pred = X_valid_reg*beta_ml + Sigma_uk*(Sigma_kk\(Y-X_reg*beta_ml));
wrong_ML = (norm(Y_valid-Y_pred)^2)/100;

var_kriging_valid = Sigma_uu - Sigma_uk*(Sigma_kk\Sigma_ku) + ... 
    (X_valid_reg' - X_reg'*(Sigma_kk\Sigma_ku))'...
    *(((X_reg'*(Sigma_kk\X_reg)))...
    \(X_valid_reg' - X_reg'*(Sigma_kk\Sigma_ku)));


X_kriging_reconstruction = [X(:,1:2); X_grid_U(:,1:2)];
D = distance_matrix(X_kriging_reconstruction);
Sigma = matern_covariance(D,par_ml(1),par_ml(2),par_ml(3));
Sigma = Sigma + par_ml(4)*eye(length(X_kriging_reconstruction));
Sigma_kk = Sigma(1:500,1:500);
Sigma_uu = Sigma(501:end,501:end);
Sigma_ku = Sigma(1:500,501:end);
Sigma_uk = Sigma(501:end,1:500);

%Field of USA_ml
Y_grid2 = nan(sz_grid);
Y_grid2(I_land) = X_grid_U*beta_ml+Sigma_uk*(Sigma_kk\(Y-X_reg*beta_ml));

var_kriging = Sigma_uu - Sigma_uk*(Sigma_kk\Sigma_ku) + ... 
    (X_grid_U' - X_reg'*(Sigma_kk\Sigma_ku))'...
    *(((X_reg'*(Sigma_kk\X_reg)))...
    \(X_grid_U' - X_reg'*(Sigma_kk\Sigma_ku)));

Y_grid2_low = Y_grid2;
Y_grid2_hi = Y_grid2;
Y_grid2_low(I_land) = X_grid_U*beta_ml+Sigma_uk*(Sigma_kk\(Y-X_reg*beta_ml)) - 1.96*sqrt(diag(var_kriging));
Y_grid2_hi(I_land) = X_grid_U*beta_ml+Sigma_uk*(Sigma_kk\(Y-X_reg*beta_ml)) + 1.96*sqrt(diag(var_kriging));



figure
imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], abs(Y_grid2-Y_grid), ...
  'alphadata', reshape(I_land,sz_grid))
axis xy tight; colorbar
title('Diff')



figure
subplot(221)
imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid2, ...
  'alphadata', reshape(I_land,sz_grid))
hold on
scatter(X(:,1), X(:,2), 20, Y, ...
  'filled','markeredgecolor','k')
scatter(X_valid(:,1), X_valid(:,2), 20, Y_valid, ...
  'filled','markeredgecolor','r')
caxis([-30 30]);
hold off; axis xy tight; colorbar
title('Kriging predictions')

subplot(222)
imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid2, ...
  'alphadata', reshape(I_land,sz_grid))
caxis([-30 30]);
hold off; axis xy tight; colorbar
title('Kriging reconstruction')

subplot(223)
imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid_low, ...
  'alphadata', reshape(I_land,sz_grid))
caxis([-30 30]);

axis xy tight; colorbar
title('Kriging low')

subplot(224)
imagesc([min(X_grid(:,1)) max(X_grid(:,1))], ...
  [max(X_grid(:,2)) min(X_grid(:,2))], Y_grid_hi, ...
  'alphadata', reshape(I_land,sz_grid))
caxis([-30 30]);

axis xy tight; colorbar
title('Kriging hi')