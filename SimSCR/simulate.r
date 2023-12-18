require(nimble)

simSCRFuncs <- nimbleFunction(
	setup = function(traps, N, lambda, sigma, beta, studyPeriod, buffer = 3){
    sigma2 <- sigma^2
    traps <- as.matrix(traps, nrow = )
    if(class(traps)[1] == "numeric") {
      traps <- matrix(traps*1.0, nrow = 1, ncol = 2)
    }else {
      traps <- as.matrix(traps*1.0)
    }
    J <- nrow(traps)
    xlim <- c(min(traps[,1]), max(traps[,1])) + c(-buffer, buffer)
    ylim <- c(min(traps[,2]), max(traps[,2])) + c(-buffer, buffer)
    S <- matrix(0, nrow = N, ncol = 2)
    studyPeriod <- studyPeriod*1.0
    paths <- array(0, c(N, 100, 3))
  },
	run = function(){},
	methods = list(
	## We will cache the key values but make it possible to update them here:
  updateStudyPeriod = function(T = double(0)){
    studyPeriod <<- T
  },
  updatePopSize = function(popSize = integer()){
    N <<- popSize
    setSize(S, c(N, 2))
  },
  updateTraps = function(newTraps = double(2)){
    J <<- dim(newTraps)[1]
    traps <<- newTraps
    xlim <<- c(min(traps[,1]), max(traps[,1])) + c(-buffer, buffer)
    ylim <<- c(min(traps[,2]), max(traps[,2])) + c(-buffer, buffer)    
  },
  updateDetFunc = function(newLambda = double(), newSigma = double(), newBeta = double(0, default = Inf)){
    sigma <<- newSigma
    sigma2 <<- sigma^2
    lambda <<- newLambda
    beta <<- newBeta
  },
  halfnormal = function(s = double(1)){
    d2 <- (traps[,1] - s[1])^2 + (traps[,2] - s[2])^2
    ans <- exp(-0.5*d2/sigma2)
    returnType(double(1))
    return(ans)
  },
	simOUPath = function(s = double(1), deltaT = double(0, default = 0.1)){
      n <- ceiling(studyPeriod/deltaT)
			mu <- matrix(value = 0, nrow = n, ncol = 3)
			rho <- exp(-beta*deltaT)
			sigmasq2 = sigma*sigma*2.0
			g <- (1-exp(-2*beta*deltaT))		
			stepSD <- sqrt(g)*sigma

			sigmasq2 = sigma*sigma*2.0
			mu[1:n,3] <- (1:n)*deltaT
			
			mu[1, 1:2] <- rnorm(2, s[1:2], sd = sigma)
			
			mu[2:n,1] <- rnorm(n-1, s[1], sd = stepSD)
			mu[2:n,2] <- rnorm(n-1, s[2], sd = stepSD)
			for(i in 2:n) mu[i,1:2] <- mu[i,1:2] + (mu[i-1, 1:2] - s[1:2])*rho	## Simulate OU process.

			returnType(double(2))
			return(mu)
	},
  simAC = function(){
    S[,1] <<- runif(N, xlim[1], xlim[2])
    S[,2] <<- runif(N, ylim[1], ylim[2])
  },
  simSCRVanillaCounts = function(newAC = logical(0, default = TRUE)){
    if(newAC) simAC()
    counts <- matrix(0, nrow = N, ncol = J)
    for(i in 1:N){
      H <- studyPeriod*lambda*halfnormal(S[i,1:2])
      for(j in 1:J){
        counts[i,j] <- rpois(1, H[j])
      }
    }
    returnType(double(2)) 
    return(counts)
  },
  simSCRVanillaOccCounts = function(nOcc = integer()){
    simAC()
    counts <- array(0, c(N, J, nOcc))
    for(k in 1:nOcc){
      counts[1:N, 1:J, k] <- simSCRVanillaCounts(FALSE)
    }
    returnType(double(3))
    return(counts)
  },
  pointInsideCircle = function(pts = double(2), j = integer(), radius = double()){
    r2 <- radius^2
    d2 <- (traps[j,1] - pts[,1])^2 + (traps[j,2] - pts[,2])^2
    ptMatch <- (d2 < r2)
    returnType(logical(1))
    return(ptMatch)
  },
  ## For this simulation, I'll assume that the data generating process is that the
  ## Animal is within the visible distance of the camera to be seen. Then it is observed
  ## with rate lambda, if the location is within the time step of the
  simContinuousMovementSCR = function(radius = double(0, default = 0.1), deltaT = double(0, default = 0.1)){
    ## Adjust lambda for the radius area.
    lamAdj <- pi*radius^2*lambda
    simAC()
    nObs <- 0
    nPhotos <- matrix(rpois(J*N, lamAdj*studyPeriod), nrow = N, ncol = J)
    times <- matrix(0, nrow = sum(nPhotos), ncol = 3)
    n <- ceiling(studyPeriod/deltaT)
    setSize(paths, c(N,n, 3))
    for( i in 1:N){
      path <- simOUPath(S[i,], deltaT)
      paths[i,,] <<- path
      for( j in 1:J ){
        timesj <- runif(nPhotos[i,j], 0, studyPeriod)
        ptsInTrap <- pointInsideCircle(path, j, radius)
        if(sum(ptsInTrap) > 0){
          tPresent <- path[ptsInTrap,3]
          for(k in seq_along(tPresent)) {
            keep <- (timesj >= tPresent[k]) & (timesj < (tPresent[k]+deltaT))
            if(any(keep)) {
              nKeep <- sum(keep)
              nObs <- nKeep + nObs
              times[(nObs-nKeep+1):nObs, 1] <- i  ## Save animal ID
              times[(nObs-nKeep+1):nObs, 2] <- j  ## Save trap
              times[(nObs-nKeep+1):nObs, 3] <- timesj[keep]  ## Save time
            }
          }
        }
      }
    }
    returnType(double(2))
    if(nObs == 0) return(matrix(0))
    return(times[1:nObs,])
  },
  getPaths = function(ID = integer()){
    returnType(double(2))
    return(paths[ID,,])
  }
))

## Build the simulation Object:
traps <- expand.grid(x = -5:5, y = -5:5)
N <- 50
lambda <- 1.5
sigma <- 2
beta <- 0.5
studyPeriod <- 20
buffer <- 5
simSCR_R <- simSCRFuncs(traps, N, lambda, sigma, beta, studyPeriod, buffer)
## Compile Object into C++
simSCR <- compileNimble(simSCR_R)

## Simulate some data:
counts <- simSCR$simSCRVanillaCounts()
## Simulate data on multiple occs. (StudyPeriod is per occ)
counts <- simSCR$simSCRVanillaOccCounts(3)

## Added an animal movement simulation. Keen to discuss this with the group...
## Playing around with it there seem to be some very silly looking results...
lambdaAdj <- 10*studyPeriod*pi^2*0.2
simSCR$updateDetFunc(newLambda = 10, newSigma = sigma, newBeta = 1)  ## Lambda should be pretty big here per unit.
times <- simSCR$simContinuousMovementSCR(radius = 0.25, deltaT = 0.01)
id <- as.numeric(names(which.max(table(times[,1]))))
sum(times[,1] == id)
plot(traps, pch = 16)
points(traps[times[times[,1] == id,2],], pch = 16, col = 'red')
path <- simSCR$getPaths(id)
lines(path[,1:2], col = 'blue')