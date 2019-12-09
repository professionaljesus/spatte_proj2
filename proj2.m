clear
close all
load("HA2_Brazil.mat")

%extract observations, E and covariates
Y = Insurance(:,2);
E = Insurance(:,1);
B = Insurance(:,4:end);

%observation matrix for all locations
A = speye(length(Y));
%find missing (nan) observations
I = ~isnan(Y);

%we need a global variable for x_mode to reuse it
%between optimisation calls
global x_mode;
x_mode = [];
%attempt to estimate parameters (optim in optim...)
%subset to only observed points here
par = fminsearch( @(x) gmrf_negloglike_NG_skeleton(x, Y(I), A(I,:), ...
B(I,:), G, E(I)), [log(100000) log(1e-5)]);
%conditional mean is now given be the mode
E_xy = x_mode;

%use the taylor  to compute posterior precision
%you need to reuse some of the code from GMRF_negloglike_NG
%to create inputs for this function call
[~, ~, Q_xy] = gmrf_taylor_skeleton(E_xy,Y(I), [A(I,:) A(I,:) B(I,:)], blkdiag(par(1)*G,par(2)*eye(length(G)),1e-6*speye(size(B(I,:),2))));
 

