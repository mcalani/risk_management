% V0 = 0.02;
% theta= 0.02;
% kappa= 0.1;
% rho = 0;
% nu = 0.2;
% T=10;

call_bsm_cf(100, 0.20, 0.02, 1, 100)



function y = chfun_norm(s0, v, r, t, w)
% Characteristic function of BSM.
% y = chfun_norm(s0, v, r, t, w)
% Inputs:
% s0: stock price
% v: volatility
% r: risk-free rate
% t: time to maturity
% w: points at which to evaluate the function
mean =log(s0)+ (r-v^2/2)*t; % mean
var = v^2*t; % variance
y = exp((i.*w*mean)-(w.*w*var*.5)); % characteristic function of log (St) evaluated at points w
end

function y = call_bsm_cf(s0, v, r, t, k)
% BSM call value calculated using formulas 2.5 to 2.7
% y = call_bsm_cf(s0, k, v, r, t, w )
% Inputs:
% s0: stock price
% v: volatility
% r: risk-free rate
% t: time to maturity
% k: option strike
% chfun_norm: Black-Scholes characteristic function
% 1st step: calculate pi1 and pi2
% Inner integral 1
int1 = @(w,s0,v,r,t,k) real(exp(-i.*w*log(k)).*chfun_norm(s0,v,r,t,w-i)./(i*w.*chfun_norm(s0, v, r, t, -i)));
int1 = integral(@(w)int1(w,s0,v,r,t,k),0,100); %numerical integration
pi1 = int1/pi+0.5;
% Inner integral 2
int2 = @(w,s0,v,r,t,k) real(exp(-i.*w*log(k)).*chfun_norm(s0, v, r, t, w)./(i*w));
int2 = integral(@(w)int2(w,s0, v, r, t, k),0,100); %numerical integration
pi2 = int2/pi+0.5; % final pi2
% 2nd step: calculate call value
y = s0*pi1-exp(-r*t)*k*pi2;
end

