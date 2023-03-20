clearvars
clc
%Parameters

sigma = 1.50; % utility parameter
delta = 0.10; % depreciation rate
beta = 0.95; % discount factor
alpha = 0.30; % capital elasticity of output
rho = 0.80; % persistence of the shock
se = 0.12; % volatility of the shock
nbk = 1000; % number of data points in the grid
nba = 2; % number of values for the shock
crit = 1; % convergence criterion
epsi = 1e-6; % convergence parameter
% 
%Discretization of the shock
% 
p = (1+rho)/2;
PI = [p 1-p;1-p p];
se = 0.12;
ab = 0;
a1 = exp(-se*se/(1-rho*rho));
a2 = exp(se*se/(1-rho*rho));
A = [a1 a2];
% 
%Discretization of the state space
%
kmin = 0.2;
kmax = 6;
kgrid = linspace(kmin,kmax,nbk)';
c = zeros(nbk,nba);
kp0 = kgrid; % initial guess on k(t+1)
util = zeros(nbk,nba);
v = zeros(nbk,nba);
Tv = zeros(nbk,nba);
Ev = zeros(nbk,nba); % expected value function
dr0 = repmat([1:nbk]',1,nba); % initial guess
dr = zeros(nbk,nba); % decision rule (will contain indices)
% 
%Main loop
%
% size(kgrid)

iter=1;
while crit>epsi
    for i=1:nbk
        for j=1:nba
            c = A(j)*kgrid(i)^alpha+(1-delta)*kgrid(i)-kgrid;
            neg = find(c<=0);
            c(neg) = NaN;
            util(:,j) = (c.^(1-sigma)-1)/(1-sigma);
            util(neg,j) = -inf;
            Ev(:,j) = v*PI(j,:)';
        end
        [v1,dr(i,:)]= max(util+beta*Ev);
    end
    % 
    %    decision rules
    %
    kp = kgrid(dr);
    Q = sparse(nbk*nba,nbk*nba);
    for j=1:nba
        c = A(j)*kgrid.^alpha+(1-delta)*kgrid-kp(:,j);
        % 
        %update the value
        %
        util(:,j)= (c.^(1-sigma)-1)/(1-sigma);
        Q0 = sparse(nbk,nbk);
        for i=1:nbk
            Q0(i,dr(i,j)) = 1;
        end
        Q((j-1)*nbk+1:j*nbk,:) = kron(PI(j,:),Q0);
    end
    Tv = (speye(nbk*nba)-beta*Q)\util(:);
    crit= max(max(abs(kp-kp0)))
    v = reshape(Tv,nbk,nba);
    kp0 = kp;
    iter= iter+1;
end
    c=zeros(nbk,nba);
for j=1:nba;
    c(:,j) = A(j)*kgrid.^alpha+(1-delta)*kgrid-kp(:,j);
end


%% PLOTS
% spy(Q)
size((speye(nbk*nba)-beta*Q)\util(:))
size(inv(speye(nbk*nba)-beta*Q))
size(util)

aux1=[2 5 6 3; 3 4 2 3 ; 2 3 4 3; 4 6 7 6]
aux2=[ 3 4 ; 5 6 ]

aux1\aux2

figure
subplot(1,3,1)
plot(kgrid,v)
title('Value Function','FontSize',14, 'interpreter','latex')
subplot(1,3,2)
plot(kgrid,kp)
title('Policy Function k','FontSize',14, 'interpreter','latex')
subplot(1,3,3)
plot(kgrid,c)
title('Policy Function c','FontSize',14, 'interpreter','latex')








