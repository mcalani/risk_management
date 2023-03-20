function ret = priceByMonteCarlo (...)
paths = generatePricePaths (...);
payoffs = computeOptionPayoffs (...);
ret = mean ( payoffs )* exp(-r*T);
end