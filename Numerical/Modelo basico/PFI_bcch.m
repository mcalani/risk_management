clearvars
clc
%Parameters

sigma = 1.50; % utility parameter
delta = 0.10; % depreciation rate
beta = 0.95; % discount factor
alpha = 0.30; % capital elasticity of output
nbk = 1000; % number of data points in the grid
crit = 1; % convergence criterion
epsi = 1e-6; % convergence parameter
ks = ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));
dev = 0.9; % maximal deviation from steady state
kmin = (1-dev)*ks; % lower bound on the grid
kmax = (1+dev)*ks; % upper bound on the grid
kgrid = linspace(kmin,kmax,nbk)'; % builds the grid
v = zeros(nbk,1); % value function
kp0 = kgrid; % initial guess on k(t+1)
dr = zeros(nbk,1); % decision rule (will contain indices)


% 
%Main loop
%
while crit>epsi
    for i=1:nbk
        % 
        %        compute indexes for which consumption is positive
        %
        imax = min(floor((kgrid(i)^alpha+(1-delta)*kgrid(i)-kmin)/dev)+1,nbk);
        % 
        %        consumption and utility
        % 
        c    = kgrid(i)^alpha+(1-delta)*kgrid(i)-kgrid(1:imax);
        util = (c.^(1-sigma)-1)/(1-sigma);
        % 
        %        find new policy rule
        %
        [v1,dr(i)]= max(util+beta*v(1:imax));
    end
    % 
    %decision rules
    %
    kp = kgrid(dr);
    c = kgrid.^alpha+(1-delta)*kgrid-kp;
    % 
    %update the value
    %
    util= (c.^(1-sigma)-1)/(1-sigma);
    Q = sparse(nbk,nbk);
    for i=1:nbk
        Q(i,dr(i)) = 1;
    end
    Tv = (speye(nbk)-beta*Q)\util;
    crit= max(abs(kp-kp0))
    v = Tv;
    kp0 = kp;
end


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
%print -dpng results.png
