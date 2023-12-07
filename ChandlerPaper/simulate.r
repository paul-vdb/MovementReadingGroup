require(nimble)

simfuncs <- nimbleFunction(
	setup = function(){},
	run = function(){},
	methods = list(
	## Do the same cumulative path calculate but for the Chandler et al. universe.
	calcChandlerPath = function(x = double(2), mu = double(2), deltaT = double(),
						sigma = double(), beta = double(), s = double(1)){
			
			N <- dim(mu)[1]
			J <- dim(x)[1]
			path <- numeric(value = 0, length = J)
			sigmasq2 = sigma*sigma*2.0
			for( i in 1:N ){
					d2 <- (x[1:J, 1] - mu[i,1])^2 + (x[1:J, 2] - mu[i,2])^2
					path <- path + exp(-d2/sigmasq2)*deltaT
			}
			returnType(double(1))
			return(path)
	},
	## Simulate a movement path.
	simPath = function(n = double(), s = double(1), sigma = double(), beta = double(), deltaT = double()){

			mu <- matrix(value = 0, nrow = n+1, ncol = 3)
			rho <- exp(-beta*deltaT)
			sigmasq2 = sigma*sigma*2.0
			g <- (1-exp(-2*beta*deltaT))		
			stepSD <- sqrt(g)*sigma

			sigmasq2 = sigma*sigma*2.0
			mu[1:(n+1),3] <- (0:n)*deltaT
			
			mu[1, 1:2] <- rnorm(2, s[1:2], sd = sigma)
			
			mu[2:(n+1),1] <- rnorm(n, s[1], sd = stepSD)
			mu[2:(n+1),2] <- rnorm(n, s[2], sd = stepSD)
			for(i in 2:(n+1)) mu[i,1:2] <- mu[i,1:2] + (mu[i-1, 1:2] - s[1:2])*rho	## Simulate OU process.

			returnType(double(2))
			return(mu)
	}	
))