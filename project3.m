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
        up = i;
    end
    
    down = i + 1;
    if mod(i, 87) == 0
        down = i;
    end
    
    left = i - 87;
    if i < 88
        left = i;
    end
    
    right = i + 87;
    if i > 8787
        right = i;
    end
    
    G(i,up) = G(i,up) - 1;
    G(i,down) = G(i,down) - 1;
    G(i,left) = G(i,left) - 1;
    G(i,right) = G(i,right) - 1;

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
