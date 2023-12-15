## Here I want to compare a simulated detection path from Chandler et al. 2021 
## and from the new method. All my helper functions are in a nimbleFunction method
## and need to be called and compiled...

library(nimble)
library(tidyverse)
source("simulate.R")

## Set up the function.
setup <- simfuncs()
simulate <- compileNimble(setup)

## Choose some parameter values from Chandler:
rho <- 0.952
beta <- -log(rho)
sigma <- exp(6.49)
lambda <- 0.952 #per hour
sigma_sm <- 26.78   ## From paper sigma_det
# sigma_sm <- sqrt(sigma^2*(1-rho^2))   ## What Paul would do (but without telemetry data).
Time <- 24 # One day simulated path.
s <- c(0,0)
mu <- simulate$simPath(n = Time, s = s, sigma = sigma, beta = beta, deltaT = 1)
plot(mu[,1:2], type = 'l')

grd <- as.matrix(expand.grid(x = seq(-1500, 1500, length = 100), y = seq(-1500,1500, length = 100)))

H.chand <- simulate$calcChandlerPath(x = grd, mu = mu[,1:2], deltaT = 1,
						sigma = sigma_sm, beta = beta, s = s)
							
dat <- data.frame(grd, H = H.chand)							

## Visualize what happens
## This is a plot of the sum of the daily intensity from a half-normal at each mask point conditional 
## on the animal's location for that hour.
ggplot(data = dat, aes(x=x, y=y)) + 
	geom_tile(aes(fill = H)) + 
	geom_path(data = data.frame(mu), aes(x=X1, y=X2), col = 'red') + 
  geom_point(data = data.frame(mu), aes(x=X1, y=X2), col = 'black')