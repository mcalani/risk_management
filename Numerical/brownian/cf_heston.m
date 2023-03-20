function y = cf_heston(u,lnS,T,r,d,V0,theta,kappa,omega,rho)

alpfa= -.5*(u.*u + u*li);
beta = kappa - rho*omega*u*li;
omega2 = omega*omega;
gamma = .5 * omega2;

D = sqrt(beta.* beta - 4.0*alpha .* gamma);
bD = beta - D;
eDt = exp(-D * T);
G = bD ./ (beta+D);
B = (bD ./ omega2) .* ((1.0 - eDt)./(1.0 - G -* eDt));
psi = (G.* eDt -1.0) ./ (G-1.0);
A = (kappa*theta)/(omega2) * (bD*T - 2.0*log(psi));

y= A + B*V0 + li*u*(lnS + (r-d) * T);

end