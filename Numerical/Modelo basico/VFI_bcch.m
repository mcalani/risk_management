clearvars
clc
%Parameters
sigma = 2;
beta= 0.95; 
delta= 0.1;
delta2= 0.05;
alpha= 0.33;
kstar=(alpha/((1/beta)-(1-delta)))^(1/(1-alpha));
kmin = 0.25*kstar; % minimum
kmax = 1.75*kstar; % maximum
kgrid = 120;
N=kgrid;
kmat=linspace(kmin,kmax,kgrid);
kmat=kmat';
v0 = zeros(kgrid,1);
v1=zeros(kgrid,1);
c = zeros(kgrid,1);
maxit = 400; % maximum number of iterations
val=zeros(kgrid,1);
tol = 100000; 
k11=zeros(kgrid,1);
iter=1;

while tol >1e-5
for i = 1:1:N
    k0=kmat(i,1);
    for j=1:1:N
        c(j,1) = k0^alpha -kmat(j,1) + (1-delta)*k0; % consumption
        if c(j,1)<=0
            val(j,1) = -9999999 - 999*abs(c(j,1));
        else
            val(j,1) = (1/(1-sigma))*(c(j,1)^(1-sigma)-1) + beta*v0(j,1);
        end
    end
v1(i,1) = max(val);
Indices = find(val==v1(i,1)); 
k11(i,1) = kmat(Indices,1);
end

tol = norm(v1-v0);
v0 = v1;
if iter>maxit
        disp('Value Function Not Converged')
        break
 end
iter = iter+1 ; 
end

V=v0;
for i = 1:1:N
c(i,1)=kmat(i,1).^alpha+(1-delta).*kmat(i,1)-k11(i,1);
end


figure
subplot(1,3,1)
plot(V)
title('Value Function','FontSize',14, 'interpreter','latex')
subplot(1,3,2)
plot(k11)
title('Policy Function k','FontSize',14, 'interpreter','latex')
subplot(1,3,3)
plot(c)
title('Policy Function c','FontSize',14, 'interpreter','latex')
print -dpng results.png


