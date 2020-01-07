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
sz = size(img);
Akron = kron(Xa, speye(sz(1)*sz(2)));
%Gibb loop

A = speye(8874);

%init values
tq = [0.1 0.1 0.1];
te = 0.1*ones(size(Y));


Nim = 2000;
for i = 1:Nim
    
    
    
    
    
    Q = kron(diag(tq), G);
    p = amd(Q);
    A_p = A(:,p);
    R = chol(Q(p,p));
    
    Q_e = diag(kron(ones(length(img),1),
    E_XY =
    Q_XY =
    
    X = R \ randn(size(R,1),1);
    EX = R \ (R' \ (A_p'*Y));
    
    
    
end



%Posterior expectation/variance
% mean(beta) ?? std(beta)?????


%Significant pixel
