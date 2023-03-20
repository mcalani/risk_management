import numpy as np

# simulate price path of lognormal distribution
p = 0.5
p_0 = 100 # Initial Price
mu = 0.06
sigma = 0.4
Np = 100000 # number of paths
Nt = 252 # number of samples
dt = 1/Nt # Time step

z = np.sign(p - np.random.rand(Nt, Np))  # discrete random variable z=+/-1 with probability p=0.5
r = mu*dt + sigma*np.sqrt(dt)*z          # dynamics/model
r = np.insert(r, 0, np.zeros(Np), axis = 0) # add zero at the start of each path
p = p_0*np.exp(r.cumsum(axis = 0)) # exponentiate the cumulative sum of returns to get the prices

# Calculate mean & std of prices at terminal Time (use p[-1] to get last price in all paths)
mean_of_distribution = p[-1].mean()
standard_deviation_of_distribution = p[-1].std()

print("Mean of distribution = {:.3f}".format(mean_of_distribution))
print('Standard deviation of distribution = {:.3f}'.format(standard_deviation_of_distribution))