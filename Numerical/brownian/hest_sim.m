nPaths=1000;
nSteps=500;
T = 1;
S0 = 1;
kappa = 2;
theta = 0.04;
v0 = 0.04;
rho = 0.99;
xi=0;

hest_sims=generatePricePathsHeston(S0, v0, kappa, theta, xi, rho, T, nPaths, nSteps);

plot(hest_sims')

