% Load data
clear
close all
load fmri.mat


col_img = colstack(img);
G = 4*eye(8874,8874);

up = -1;
down = 1;
left = 0;
right = 0;
for i = 1:8874
    up = i - 1;
    if mod(i,88) == 0 || i == 1
        G(i,i) = G(i,i) - 1;
    else
        G(i,up) = G(i,up) - 1;
    end    
    down = i + 1;
    if mod(i, 87) == 0
        G(i,i) = G(i,i) - 1;
    else
        G(i,down) = G(i,down) - 1;
    end
    left = i - 87;
    if i < 88
        G(i,i) = G(i,i) - 1;
    else
        G(i,left) = G(i,left) - 1;
    end
    right = i + 87;
    if i > 8787
        G(i,i) = G(i,i) - 1;
    else
        G(i,right) = G(i,right) - 1;
    end
end
G = sparse(G);




Xa = [X(:,1:2) sum(X(:,3:end),2)];

Y = colstack(col_img);
A = kron(Xa,speye(8874,8874));

%Gibb loop

%init values
tq = [0.1 0.1 0.1];
te = 0.1*ones(8874,1);
Nim = 100;
tq_hist = zeros(3, Nim);
te_hist = zeros(8874, Nim);
beta_hist = zeros(26622, Nim);
burnin = 10;

for i = 1:Nim
    Q = kron(sparse(diag(tq)), G);
    Q_e = spdiags(kron(ones(160,1),te),0,1419840,1419840);
    Q_xy = Q+A'*Q_e*A;
    
    p = amd(Q_xy);
    R = chol(Q_xy(p,p));
    A_p = A(:,p);
    Y_p = Y(p);
    
%     EX = R \ (A_p'*Q_e*Y_p);
    EX = R\((A_p'*Q_e*Y)'/R)';
    beta = EX + R\randn(size(R,1),1);
    beta(p) = beta;
    beta_hist(:, i) = beta;
    
    % tq1
    N = length(G);
    shape = N/2 + 1;
    scale = 2/(beta(1:8874)'*G*beta(1:8874));
    tq_1 = gamrnd(shape, scale);
    
    %tq2
    scale = 2/(beta(8875:17748)'*G*beta(8875:17748));
    tq_2 = gamrnd(shape, scale);
    
    %tq3
    scale = 2/(beta(17749:end)'*G*beta(17749:end));
    tq_3 = gamrnd(shape, scale);
    tq = [tq_1 tq_2 tq_3];
    tq_hist(:,i) = tq;
    
%     te / test
    N = 160;
    shape = N/2 + 1;

    e_sample = Y-A*beta;
    for s = 1:8874
        scale = 2/(e_sample(s:8874:end)'*e_sample(s:8874:end));
        te(s,1) = gamrnd(shape, scale);
    end    
    te_hist(:,i) = te;

    i
end
%%
beta_mean = mean(beta_hist(:, burnin:end),2);
beta_recon = reshape(beta_mean, [87, 102 ,3]);

te_mean = mean(te_hist(:, burnin:end),2);
tq_mean = mean(tq_hist(:, burnin:end),2);
Q_recon = kron(sparse(diag(tq_mean)), G);
Q_e_recon = spdiags(kron(ones(160,1),te_mean),0,1419840,1419840);
Q_xy_recon = Q_recon+A'*Q_e_recon*A;

Y_recon = A*beta_mean;
img_recon = reshape(Y_recon, [87, 102, 160]);

diff = img - img_recon;

%Variance boi
p = amd(Q_xy_recon);
R = chol(Q_xy_recon(p,p));


beta_variance = zeros(length(R), 100);
for i = 1:100
    X = R\randn(size(R,1),1);
    beta_variance(:,i) = X(p);
end
beta_variance = mean(beta_variance,2);
beta3_mean = beta_mean(end - 8873:end,1);

beta3_variance = beta_variance(end - 8873:end, 1);
beta3_variance_recon = reshape(beta3_variance, [87, 102]);


sig_img = ones(8874,1);
for jj = 1:8874
    if beta3_mean(jj,1) - 4.8272*sqrt(beta3_variance(jj,1))/sqrt(90) < 0 && ...
            beta3_mean(jj,1) + 4.8272*sqrt(beta3_variance(jj,1))/sqrt(90) > 0
        sig_img(jj) = 0;
    end
end

sig_img = reshape(sig_img, [87, 102]);

%%

%Figures

figure
imagesc(beta_recon(:,:,3))
title("Beta_3")

figure
for t = 1:160
    subplot(131)
    imagesc(img(:,:,t))
    subplot(132)
    imagesc(img_recon(:,:,t))
    subplot(133)
    imagesc(diff(:,:,t))
    colorbar
    drawnow
    pause(0.1)
end
