function [ prices, variances ] = generatePricePathsHeston( ...
    S0, v0, ...
    kappa, theta, xi, rho, ...
    T, nPaths, nSteps)
    %GENERATEPRICEPATHSHESTON Generate price paths according to the
    % Heston model
    prices = zeros( nPaths, nSteps );
    variances = zeros( nPaths, nSteps );
    currS = S0;
    currv = v0;
    dt = T/nSteps;
    for i=1:nSteps
    %epsilon = mvnrnd( [0 0],[1 rho; rho 1], nPaths );
    epsilon = randnMultivariate( [1 rho; rho 1], nPaths );
    dW1 = epsilon(1,:)*sqrt(dt);
    dW2 = epsilon(2,:)*sqrt(dt);
    currS = currS + sqrt( currv).* currS .* dW1';
    currv = currv + kappa*(theta - currv)*dt + xi*sqrt( currv).* dW2';
    currv = abs( currv ); % Forcibly prevent negative variances
    prices( :, i) = currS;
    variances( :, i) = currv;
end