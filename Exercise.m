
N1 = [0 1 0; 1 0 1; 0 1 0];
N2 = [1 1 1; 1 0 1; 1 1 1];
beta = 2;
alpha = 0.1;
z0 = randn(120,120,2);
[z,Mz,Mf] = mrf_sim(z0,N2,alpha,beta,100);

for iter=1:10
    z = mrf_sim(z,N1,alpha,beta,100);
    image(rgbimage(z))
    drawnow
end


% [z,Mz,Mf,Mzf]=mrf_sim(z0,N,alpha,beta,iter,gt)
%   N = [0 1 0;1 0 1;0 1 0];
%   z = mrf_sim(zeros(100,120,3),N,log([0.2 0.3 0.5]),0.5,100);
%   for iter=1:100
%     z = mrf_sim(z,N,log([0.2 0.2 0.6]),[0.9 0.9 0],1);
%     image(rgbimage(z))
%     drawnow
%   end
