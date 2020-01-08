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
        G(i,i) = G(i,i) - 0;
    else
        G(i,up) = G(i,up) - 1;
    end    
    down = i + 1;
    if mod(i, 87) == 0
        G(i,i) = G(i,i) - 0;
    else
        G(i,down) = G(i,down) - 1;
    end
    left = i - 87;
    if i < 88
        G(i,i) = G(i,i) - 0;
    else
        G(i,left) = G(i,left) - 1;
    end
    right = i + 87;
    if i > 8787
        G(i,i) = G(i,i) - 0;
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
tau_hist = zeros(Nim,3);
epsilon_hist = zeros(Nim,8874);
beta_hist = zeros(Nim,26622);

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
    x_samp = EX + R\randn(size(R,1),1);
    x_samp(p) = x_samp;
    
    % tq1
    N = length(G);
    shape = N/2 + 1;
    scale = 2/(x_samp(1:8874)'*G*x_samp(1:8874));
    tq_1 = gamrnd(shape, scale);
    
    %tq2
    scale = 2/(x_samp(8875:17748)'*G*x_samp(8875:17748));
    tq_2 = gamrnd(shape, scale);
    
    %tq3
    scale = 2/(x_samp(17749:end)'*G*x_samp(17749:end));
    tq_3 = gamrnd(shape, scale);
    tq = [tq_1 tq_2 tq_3];
    tau_hist(i,:) = tq';
    
%     te / test
    N = length(Q_e);
    shape = N/2 + 1;

    e_sample = Y-A*x_samp;
    for s = 1:8874
        scale = 2/(e_sample(s:8874:end)'*e_sample(s:8874:end));
        te(s,1) = gamrnd(shape, scale);
        epsilon_hist(i,:) = te';

    end    
    i
end

TT = reshape(x_samp, [87, 102,3]);


%Posterior expectation/variance
% mean(beta) ?? std(beta)?????


%Significant pixel