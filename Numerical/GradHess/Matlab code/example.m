% Code to
% a) simulate linear model, y_i = 0.2 + 0.2 x_1i -0.1 x_2i + 0.9 x_3i + err_i
%       x_1i and x_s1 come from N(0,1)
%       err_i ~ N(0,sd^2)
% b) estimate it by OLS
% c) estimate ot by ML
% d) estimate different ML standard errors
% e) Test multiple restriction
%
% This code requires the following Functions:
% OLSest, nll_lin, HessMp, gradp
clc
clear all

%% a) simulate linear model

T   = 100;                     % set sample size
b0  = [0.2; 0.2; -0.1; 0.8; 0];      % set parameter values
sd  = 1.0;                      % set error standard deviation
x   = [ones(T,1) randn(T,4)];   % simulate X matrix
err = randn(T,1)*sd;            % simulate error terms
y = x*b0 + err;                 % calculate y_i s

%% b) estimate it by OLS
 b=(x'*x)\(x'*y) ;
 y_hat    = x*b;
 err      = y - y_hat;
 sigma    = (err'*err)/T;
 var_cov  = inv(x' *x) *  sigma;
 bse      = sqrt(diag(var_cov));
 
%[b,bse,res,n,rss,r2] = OLSest(y,x,1);   % OLS estimation

%% c) estimate ot by ML
datamat = [y x];                        % define data matrix for use in nll_lin
theta0 = [mean(y); zeros(size(b0,1)-1,1); std(y)];       % this sets the initial parameter vector
options = optimset; % sets optimisation options to default
[thetaopt] = fminunc(@nll_lin,theta0,options,datamat,1);

%% d) estimate different ML standard errors

H = HessMp(@nll_lin,thetaopt,datamat,1); % this returns the negative of the Hessian
g = Gradp(@nll_lin,thetaopt,datamat,0);  % this returns a (T x size(thetaopt,1)) matrix of gradients
J = (g'*g);                             % calculates the OPG

se_H = sqrt(diag(inv(H)));
se_J = sqrt(diag(inv(J)));
se_SW = sqrt(diag(inv(H*inv(J)*H)));    % Sandwich variance covariance

disp('   Est     se(OLS)     se(H)     se(J)     se(SW)');
disp([thetaopt [bse;0] se_H se_J se_SW]);

%% e) Testing multiple restriction
% b(2) = 0, b(3) = 0

% LR test
% estimate the restricted model
theta0_r = [mean(y); 0; 0; std(y)];       % this sets the initial parameter vector (4x1) 
        % do not hand in x1 and x2
[theta_r] = fminsearch(@nll_lin,theta0_r,options,datamat(:,[1 2 5 6]),1);

L_u = -nll_lin(thetaopt,datamat,1);             % calculates unrestricted logLikelihood
L_r = -nll_lin(theta_r,datamat(:,[1 2 5 6]),1); % calculates restricted logLikelihood

LR   = 2*(L_u - L_r);
LR_p = 1-chi2cdf(LR,2);

fprintf('LR test = %6.2f; p-value = %6.4f \n', LR, LR_p);


% Wald Test
% Specify restriction matrix R
R = [0 1 0 0 0 0; 0 0 1 0 0 0]; % (2x6) restriction matrix 
b = [0;0];                  % (2x1) constraints

V = inv(H);

rest = (R*thetaopt - b);

W = rest'*inv(R*V*R')*rest;
W_p = 1-chi2cdf(W,2);

fprintf('Wald test = %6.2f; p-value = %6.4f \n', W, W_p);

% LM test

% construct restricted (but full) parameter vector
theta_r_full = [theta_r(1); 0; 0; theta_r(2:end)];

G = Gradp(@nll_lin,theta_r_full,datamat,1);  % this returns a (1 x size(thetaopt,1)) vector of gradients

H_r = HessMp(@nll_lin,theta_r_full,datamat,1); % this returns the negative of the Hessian at theta_r_full
V_r = inv(H_r);                                % calculates V(theta_r_full)

LM = G*V_r*G';
LM_p = 1-chi2cdf(LM,2);

fprintf('LM test = %6.2f; p-value = %6.4f \n', LM, LM_p);