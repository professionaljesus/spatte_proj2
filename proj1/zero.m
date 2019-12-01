load T_lund.mat


t = T_lund(:,1);  Y = T_lund(:,2);  n = length(Y);
X = [ones(n,1) sin(2*pi*t/365) cos(2*pi*t/365)];
beta = regress(Y, X);eta = Y-X*beta;plot(t, Y, 'b', t, X*beta, 'r');datetick
alpha = regress(eta(2:end),eta(1:end - 1));
v = eta(2:end) - eta(1:end - 1)*alpha;
sigma = std(v);

t0 = 500;
esty = zeros(25,1);
vary = [sigma; zeros(24,1)];
for k = 1:25
    esty(k,1) = X(t0 + k - 1,:)*beta + (alpha^k)*eta(t0);
    if k > 1
        vary(k,1) = vary(k - 1,1) + alpha^(2*(k - 1))*sigma;
    end
end
 plot([(esty + 1.96*vary),(esty - 1.96*vary),esty,Y(t0:t0+24)])