% Load data
clear
close all
load fmri.mat

%Stack pixels
col_img = colstack(img);

%Create Neigourborhood with Dirchlet edges
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

%Create Model unchaning parameters
Xa = [X(:,1:2) sum(X(:,3:end),2)];
Y = colstack(col_img);
A = kron(Xa,speye(8874,8874));

%Gibb loop
%init values
tq = [0.1 0.1 0.1];
te = 0.1*ones(8874,1);
Nim = 200;
tq_hist = zeros(3, Nim);
te_hist = zeros(8874, Nim);
beta_hist = zeros(26622, Nim);
burnin = 10;

%Loop
for i = 1:Nim
    
    %Beta sampling
    
    %Create recounstruction of GMRF
    Q = kron(sparse(diag(tq)), G);
    Q_e = spdiags(kron(ones(160,1),te),0,1419840,1419840);
    Q_xy = Q+A'*Q_e*A;
    
    p = amd(Q_xy);
    R = chol(Q_xy(p,p));
    A_p = A(:,p);
    Y_p = Y(p);
    EX = R\((A_p'*Q_e*Y)'/R)';
    beta = EX + R\randn(size(R,1),1);
    beta(p) = beta;
    beta_hist(:, i) = beta;
    
    %Tau_Q
    
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
    
    %Tau_Epsilon
    N = 160;
    shape = N/2 + 1;
    e_sample = Y-A*beta;
    
    %Pixelvies
    for s = 1:8874
        scale = 2/(e_sample(s:8874:end)'*e_sample(s:8874:end));
        te(s,1) = gamrnd(shape, scale);
    end    
    te_hist(:,i) = te;
end
%%

%From the Gibbs, construct a mean sample
beta_mean = mean(beta_hist(:, burnin:end),2);
beta_recon = reshape(beta_mean, [87, 102 ,3]);
te_mean = mean(te_hist(:, burnin:end),2);
te_variance = std(te_hist');
tq_mean = mean(tq_hist(:, burnin:end),2);

%Reconstruct Beta, from parameters
Q_recon = kron(sparse(diag(tq_mean)), G);
Q_e_recon = spdiags(kron(ones(160,1),te_mean),0,1419840,1419840);
Q_xy_recon = Q_recon+A'*Q_e_recon*A;

% Variance 
p = amd(Q_xy_recon);
R = chol(Q_xy_recon(p,p));
A_p = A(:,p);
beta_variance = R\randn(size(R,1),1);
beta_expectation =  R\((A_p'*Q_e_recon*Y)'/R)';

%Reconstruction of picutre
Y_recon = A*beta_mean;
Y_diff = Y-Y_recon;
img_recon = reshape(Y_recon, [87, 102, 160]);
diff = img - img_recon;

%Study tau_e
% figure
% imagesc(reshape(te_mean,[87 102]));
% figure
% imagesc(reshape(te_variance',[87 102]));
% residual_with_time = reshape(Y_diff,[87 102 160]);
% residual_avreage = mean(residual_with_time(:,:, burnin:end),3);
% figure
% imagesc(residual_avreage);


%Calculate signification
I_significant = abs(beta_mean(17749:end)) > norminv(1-0.05/8775)*sqrt(abs(beta_variance(17749:end)));
signifanct_procent = sum(I_significant)/8774;

%%

%Figures

% figure
% imagesc(beta_recon(:,:,3))
% title("Beta_3")
% 
% figure
% for t = 1:160
%     subplot(131)
%     imagesc(img(:,:,t))
%     subplot(132)
%     imagesc(img_recon(:,:,t))
%     subplot(133)
%     imagesc(diff(:,:,t))
%     colorbar
%     drawnow
%     pause(0.1)
% end
