% Load data
clear
close all
load fmri.mat

Xa = [X(:,1:2) sum(X(:,3:end),2)];

col_img = colstack(img);
G = 4*speye(8874,8874);

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

tau = 1;
Q = tau * G;

Y = img(:);
sz = size(img);
A = kron(Xa, speye(sz(1)*sz(2)));
%Gibb loop



%Posterior expectation/variance
% mean(beta) ?? std(beta)?????


%Significant pixel
