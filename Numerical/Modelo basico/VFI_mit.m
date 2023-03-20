%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Simple VFI Model with MIT Shock
%  Auhor: S. Ramirez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
clc
%Parameters
sigma = 2;
beta= 0.95; 
delta= 0.15;
delta2= delta*(1.25); % Initial ss Option 1
delta3= delta*(0.75); % Initial ss Option 2
alpha= 0.33;
kstar=(alpha/((1/beta)-(1-delta)))^(1/(1-alpha));
kstar2=(alpha/((1/beta)-(1-delta2)))^(1/(1-alpha));
kstar3=(alpha/((1/beta)-(1-delta3)))^(1/(1-alpha));
kmin = 0.1*kstar; % minimum
kmax = 1.5*kstar; % maximum
kgrid = 3000;
N=kgrid;
kmat=linspace(kmin,kmax,kgrid);
kmat=kmat';

previos_guess=1; % if you don't have previuos guess put 0
%if previos_guess==1
%    load('Vguess.mat');
%    v0=V;
%else
v0 = zeros(kgrid,1);
%end
c = zeros(kgrid,1);
maxit = 500; % maximum number of iterations
val=zeros(kgrid,1);
v1=zeros(kgrid,1);
k11=zeros(kgrid,1);
tol = 100000; 
iter=1;
V=v0;
tic
while tol >1e-6
for i = 1:1:N
    k=kmat(i,1);
    for j=1:1:N
    c(j,1) = k(j,1)^alpha -kmat(j,1) + (1-delta)*k; % consumption
    if c(j,1)<=0
    val(j,1) = -9999999 - 999*abs(c(j,1));
    else
    val(j,1) = (1/(1-sigma))*(c(j,1)^(1-sigma)-1) + beta*V(j,1);
    end
    end
v1(i,1) = max(val);
Indices = find(val==v1(i,1)); 
k11(i,1) = kmat(Indices,1);
end
tol = norm(v1-V);
V = v1;
if iter>maxit
        disp('Value Function Not Converged')
        break
 end
iter = iter+1 ; 
end
toc
for i = 1:1:N
c(i,1)=kmat(i,1).^alpha+(1-delta).*kmat(i,1)-k11(i,1);
end
save('Vguess.mat','V');
figure
subplot(1,3,1)
%plot(V)
%hold on
plot(kmat,V,'b','LineWidth',2);
title('Value Function','FontSize',16, 'interpreter','latex')
xlabel('k','FontSize',10, 'interpreter','latex'); 
ylabel('V(k)','FontSize',10, 'interpreter','latex');
xlim([kmin kmax])
subplot(1,3,2)
%plot(k11)
%hold on
plot(kmat,k11,'b','LineWidth',2);
title('Policy Function for K','FontSize',16, 'interpreter','latex')
xlabel('k','FontSize',10, 'interpreter','latex'); 
ylabel('g(k)','FontSize',10, 'interpreter','latex');
xlim([kmin kmax])
subplot(1,3,3)
%plot(c)
%hold on
plot(kmat,c,'b','LineWidth',2);
title('Policy Function for C','FontSize',16, 'interpreter','latex')
xlabel('k','FontSize',10, 'interpreter','latex'); 
ylabel('c(k)','FontSize',10, 'interpreter','latex');
xlim([kmin kmax])
print -dpng results.png

%%
%MIT Shock

% Initial ss Option 1
i=1;
time1(i)=1;
k111(i)=kstar2;
tol = 100000; 
while tol >1e-6
i=i+1;
time1(i) = time1(i-1)+1;
k111(i) = interp1(kmat,k11,k111(i-1),'spline'); 
tol = norm(k111(i-1)-k111(i));
Period=time1(i);
end

t1=max(time1);
% Initial ss Option 2
i=1;
time2(i)=1;
k112(i)=kstar3;
tol = 100000; 
while tol >1e-6
i=i+1;
time2(i) = time2(i-1)+1;
k112(i) = interp1(kmat,k11,k112(i-1),'spline'); 
tol = norm(k112(i-1)-k112(i));
Period=time2(i);
end

t2=max(time2);
t=min(t1,t2);

% Figures
css_1=kstar2.^alpha+(1-delta2).*kstar2-kstar2;
vss_1=(1/(1-sigma))*((css_1^(1-sigma))-1); 

css_2=kstar3.^alpha+(1-delta3).*kstar3-kstar3;
vss_2=(1/(1-sigma))*((css_2^(1-sigma))-1); 

c1 = zeros(1,t);
v1 = zeros(1,t);
inv = zeros(1,t);

c1(1,1)=css_1;
v1(1,1)=vss_1; 



for i = 2:1:t1
c1(1,i)=k111(1,i-1).^alpha+(1-delta).*k111(1,i-1)-k111(1,i);
v1(1,i)=(1/(1-sigma))*(c1(1,i)^(1-sigma)-1);
inv(1,i)=k111(1,i)-k111(1,i-1);
end

c11=(c1./css_1-1)*100;
v11=((v1-vss_1)./abs(vss_1))*100;
k1111=(k111./kstar2-1)*100;

c2 = zeros(1,t);
v2 = zeros(1,t);
inv2 = zeros(1,t);


c2(1,1)=css_2;
v2(1,1)=vss_2; 



for i = 2:1:t2
c2(1,i)=k112(1,i-1).^alpha+(1-delta).*k112(1,i-1)-k112(1,i);
v2(1,i)=(1/(1-sigma))*(c2(1,i)^(1-sigma)-1);
inv2(1,i)=k112(1,i)-k112(1,i-1);
end

c22=(c2./css_2-1)*100;
v22=((v2-vss_2)./abs(vss_2))*100;
k1122=(k112./kstar3-1)*100;


figure
subplot(1,4,1)
plot(k111,'b','LineWidth',2)
hold on
plot(k112,'r','LineWidth',2)
title('Capital','FontSize',16, 'interpreter','latex')
xlim([2 t])
xlabel('Periods','FontSize',12, 'interpreter','latex')
subplot(1,4,2)
plot(c1,'b','LineWidth',2)
hold on
plot(c2,'r','LineWidth',2)
title('Consumption','FontSize',16, 'interpreter','latex')
xlim([2 t])
xlabel('Periods','FontSize',12, 'interpreter','latex')
subplot(1,4,3)
plot(inv,'b','LineWidth',2)
hold on
plot(inv2,'r','LineWidth',2)
title('Investment','FontSize',16, 'interpreter','latex')
xlim([2 t])
xlabel('Periods','FontSize',12, 'interpreter','latex')
subplot(1,4,4)
plot(v1,'b','LineWidth',2)
hold on
plot(v2,'r','LineWidth',2)
title('Utility','FontSize',16, 'interpreter','latex')
xlim([2 t])
xlabel('Periods','FontSize',12, 'interpreter','latex')
grid off
legend({'Low depreciation rate','High depreciation rate'},'FontSize',8, 'interpreter','latex') 
print -dpng mit_shock_1.png

figure
subplot(1,3,1)
plot(k1111,'b','LineWidth',2)
hold on
plot(k1122,'r','LineWidth',2)
title('Capital','FontSize',16, 'interpreter','latex')
xlim([2 t])
xlabel('Periods','FontSize',12, 'interpreter','latex')
subplot(1,3,2)
plot(c11,'b','LineWidth',2)
hold on
plot(c22,'r','LineWidth',2)
title('Consumption','FontSize',16, 'interpreter','latex')
xlim([2 t])
xlabel('Periods','FontSize',12, 'interpreter','latex')
subplot(1,3,3)
plot(v11,'b','LineWidth',2)
hold on
plot(v22,'r','LineWidth',2)
title('Utility','FontSize',16, 'interpreter','latex')
xlim([2 t])
xlabel('Periods','FontSize',12, 'interpreter','latex')
grid off
legend({'Low depreciation rate','High depreciation rate'},'FontSize',8, 'interpreter','latex') 
print -dpng mit_shock_2.png