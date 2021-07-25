#### Cesar Pardede ####
# Final Project
# Math 538
# 2020-12-01

library(plgp) # distance() - creates symmetric distance matrix
library(mvtnorm) # rmvnorm() - samples from MVN distn


f = function(x){
  # true function of x
  f = x/2*sin(3*x/4)
  f
}


kov = function(X, X1 = X, t = 1, l = 1){
  # covariance/kernel function
  # t: amplitude parameter
  # l: width parameter
  d = distance(X, X1)
  k = t^2 * exp(-d/(2*l^2))
  k
}


gaussianProcess = function(X, data, t, l, n_gp){
  # n_gp: number of GP's to generate
  # calculating parameters of GPs
  x_data = data$x
  y_data = data$y
  SXX = kov(X, t = t, l = l) # NOTE: UPPERCASE x 
  Sxx = kov(x_data, t = t, l = l) # note: lowercase x
  SXx = kov(X, x_data, t = t, l = l)
  Sxx_inv = solve(Sxx)
  # mu = SXx %*% Sxx_inv %*% y_data # mu vector with data
  mu = X + SXx %*% Sxx_inv %*% (y_data-x_data) # mu vector with data
  S = SXX - SXx %*% Sxx_inv %*% t(SXx) # sigma matrix w/ data
  
  # Gaussian Processes posterior
  Y_gp = rmvnorm(n_gp, mu, S)
  gp = list(Y = Y_gp, m = mu, S = S)
  gp
}


gp_plotter = function(X, data, title = 'Gaussian Process Posteriors', t = 1, l = 1, n_gp = 100, fp = T, dp = T){
  # plotting gp's
  x_data = data$x
  y_data = data$y
  gp = gaussianProcess(X, data, t, l, n_gp)
  Y_gp = gp$Y
  mu = gp$m
  S = gp$S
  
  plot(1, type = 'n', xlim = c(x_a, x_b), ylim = c(-10, 10), xlab = 'X', ylab = 'Y', main = title)
  for (i in 1:n_gp){
    lines(X, Y_gp[i, ], col = 8)
  }
  
  # 95% bounds
  var = diag(S)
  var[which(var < 0)] = 0
  sd = sqrt(var)
  q_lower = mu + qnorm(0.025, 0, sd)
  q_upper = mu + qnorm(0.975, 0, sd)
  lines(X, q_lower, lty = 2, col = 2, lwd = 2)
  lines(X, q_upper, lty = 2, col = 2, lwd = 2)
  
  # plotting true function
  plot_truth_and_data(data, fp, dp)
}


plot_truth_and_data = function(data, func_plot = T, data_plot = T, new_plot = F){
  x_data = data$x
  y_data = data$y
  if (new_plot){plot(1, type = 'n', xlim = c(x_a, x_b), ylim = c(-10, 10), xlab = 'X', ylab = 'Y', main = 'True Function and Data')}
  if (func_plot){lines(X, Y, col = 4)}
  if (data_plot){points(x_data, y_data, pch = 20, col = 4)}
}


gp_priors = function(n_prior, title = 'Gaussian Process Priors'){
  # Gaussian Process priors
  SXX = kov(X)
  plot(1, type = 'n', xlim = c(x_a, x_b), ylim = c(-10, 10), xlab = 'X', ylab = 'Y', main = title)
  for (i in 1:n_prior){
    mu = rep(0, ncol(SXX))
    # mu = -f(X)
    Y_gp_prior = rmvnorm(1, mu, sigma = SXX)
    lines(X, Y_gp_prior, col = 8)
  }
}


# True function
n = 200 # length of grid/discretization/quantization, and dimension of MVN
x_a = 0; x_b = 4*pi
X = matrix(seq(x_a, x_b, length = n), ncol = 1)
Y = f(X)
DATA = data.frame(x = X, y = Y)

# plotting gp priors
# Gaussian Processes are n-variate MVN distributions
gp_priors(1, title = '1 GP Prior')
gp_priors(10, '10 GP Priors')
gp_priors(100, '100 GP Priors')

# define "data"
n_data = 4 # 15 is max number of x's before numerical instability
set.seed(538)
x_data = as.matrix(sample(X, n_data)) # random data x
# x_data = as.matrix(seq(x_a, x_b, length = n_data)) # uniform data x
y_data = f(x_data)
data = data.frame(x = x_data, y = y_data)
data0 = data

# Posterior GP conditioned on data
plot_truth_and_data(data, new_plot = T, data_plot = T, func_plot = F)
gp_plotter(X, data, fp = F, dp = T)

# measure new data based on distribution of functions
x_new = 7
y_obs = f(x_new)
data_new = data.frame(x = x_new, y = y_obs)
data = rbind(data, data_new)
gp_plotter(X, data, fp = F, dp = T)

x_new = 9
y_obs = f(x_new)
data_new = data.frame(x = x_new, y = y_obs)
data = rbind(data, data_new)
gp_plotter(X, data, fp = F, dp = T)

x_new = seq(x_a, x_b, length = 9)
y_obs = f(x_new)
data_new = data.frame(x = x_new, y = y_obs)
gp_plotter(X, data_new, fp = F, dp = T)

# plot true function
plot_truth_and_data(DATA, func_plot = T, data_plot = F)

# plot comparisons of tuning parameters
data = data_new
gp_plotter(X, data, t = 1, l = 1, title = 't = 1, l = 1')
gp_plotter(X, data, t = 1/2, l = 1, title = 't = 1/2, l = 1')
gp_plotter(X, data, t = 1, l = 1/2, title = 't = 1, l = 1/2')
gp_plotter(X, data, t = 1/2, l = 1/2, title = 't = 1/2, l = 1/2')

# additional plots for slides
plot(1, type = 'n', xlim = c(x_a, x_b), ylim = c(-10, 10), xlab = 'X', ylab = 'Y', main = 'X_i from 0 to 4*pi')
for (i in 1:(length(X)/4)){
  abline(v = X[4*i], lty = 1, col = 4)
  points(X[4*i], 0, pch = 20, col = 2)
}

gp_priors(100, '100 GP Priors and Data (m = 0)')
plot_truth_and_data(data0, func_plot = F)
