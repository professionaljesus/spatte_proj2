% Load data
load fmri.mat

%size of data
sz = size(img);

%Option 1: regress onto indicator functions
beta = X\colstack(img)';
%reshape back to an image where the beta-coefs are in each "color"-layer
beta = reshape(beta', sz(1), sz(2), []);
%and treat the beta:s as the image
%Perform a PCA on the regression coefficients to find important components
[y_beta, ~, P_beta] = pca(colstack(beta));
y_beta = reshape(y_beta, sz(1), sz(2), []);

recon = colstack(beta)*X';
recon = reshape(recon, sz(1), sz(2), []);
eps = img - recon;

%%

%lets plot the pca components
figure
subplot(3,4,1)
semilogy(P_beta/sum(P_beta))
axis tight
for i=1:size(y_beta,3)
  subplot(3,4,i+1)
  imagesc(y_beta(:,:,i))
  title(i)
end

%Option 2: Compute SVD directly on the data (this is essentially the SVD
%example from lecture 9)
[y,V,P] = pca(colstack(img));
y = reshape(y, sz(1), sz(2), []);

%study the temporal components to find those with 20s periodicity
figure
subplot(3,4,1)
semilogy(P/sum(P))
axis tight
for i=1:11
  subplot(3,4,i+1)
  plot(V(:,i))
  axis tight
  title(i)
end

%alternatively using a periodogram (stationary/timeseries course)
figure
for i=1:8
  subplot(4,4,2*i-1)
  plot(V(:,i))
  axis tight
  title(i)
  subplot(4,4,2*i)
	periodogram(V(:,i))
	line([0.1 0.1], get(gca,'YLim'));
	xlabel(''); ylabel(''); title(i);
end

%and plot a few leading components
figure
subplot(3,4,1)
semilogy(P/sum(P))
axis tight
for i=1:11
  subplot(3,4,i+1)
  imagesc(y(:,:,i))
  title(i)
end
