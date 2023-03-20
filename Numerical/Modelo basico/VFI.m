clear all;
close all;
tic

global v0 beta delta alpha kmat k0 s kgrid

%set parameters
alpha= 0.33; %capital's share
beta= 0.95;
delta= 0.1; %depreciation
s=2;

tol=0.01;
maxits=300;
dif= tol+1000;
its=0;

%first solve for steady state
kstar = (alpha/ (1/beta-(1-delta)))^(1/(1-alpha)); %steady state k
cstar= kstar^(alpha)- delta*kstar;
istar= delta*kstar;
ystar=kstar^(alpha);

kmin= 0.25*kstar; %minimum
kmax= 1.75*kstar; %maximum
kgrid= 99; %grid points +1
grid= (kmax-kmin)/kgrid;%grid

kmat= kmin:grid:kmax;
kmat=kmat';

[N,n]= size(kmat);

v0= zeros(N,1);
v1= zeros(N,1);
k11=zeros(N,1);

while dif>tol && its<maxits
    for i = 1:N
        k0=kmat(i,1);
        k1=fminbnd(@valfun2, kmin,kmax);
        v1(i,1)= -valfun(k1);
        k11(i,1)=k1;
    end
    dif=norm(v1-v0);
    v0=v1;
    its= its+1;
end

toc

% This program gets the value function for a neoclassical growth model with
% no uncertainty and CRRA utility
function val=valfun2(k)

global v0 beta delta alpha kmat k0 s 

klo=max(sum(k>kmat),1);

khi=klo+1;

gg=v0(klo) + (k-kmat(klo))*(v0(khi)-v0(klo))/(kmat(khi)-kmat(klo));

c=k0^alpha - k + (1-delta)*k0; %consumption
if c<=0
    val= - 9999999 - 999*abs(c);
else
    val= (1/(1-s))*(c^(1-s)-1) + beta*gg;
end
val= -val; 
end










