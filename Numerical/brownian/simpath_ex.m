% Script to price an Asian put option using a monte-carlo approach.
S0 = [50 48 49]  ;       % Price of underlying today
mu = [0.03 0.06 0.04];     % expected return
sig = [0.05 0.1 0.08];     % expected vol.
corr = [1   0.7 0.65;
        0.8  1  0.65;
        0.7 0.6  1]; % correlation matrix
dt = 1/365;   % time steps
etime = 500;   % days to expiry
T = dt*etime; % years to expiry

nruns = 1000; % Number of simulated paths

% Generate potential future asset paths
S = AssetPathsCorrelated(S0,mu,sig,corr,dt,etime,nruns);

% Plot one set of sample paths
time = etime:-1:0;
plot(time,squeeze(S(:,130,:)),'Linewidth',2);
set(gca,'XDir','Reverse','Fontsize',20);
xlabel('Time to Expiry','Fontsize',20);
ylabel('Asset Price','Fontsize',20);
title('One Set of Simulated Asset Paths','Fontsize',20);
grid on
set(gcf,'Color','w');

%http://www.goddardconsulting.ca/option-pricing-monte-carlo-basket.html
%http://www.goddardconsulting.ca/matlab-monte-carlo-assetpaths-corr.html