library(gk)

# Function to compute pool of statistics of g-and-k distribution
gkStatistics <- function(obsData, nStats){
  S <- array(0, dim = c(nStats, 1))
  
  S[1] = quantile(obsData, 0.5)
  
  S[2] = quantile(obsData, 0.75) - quantile(obsData, 0.25)
  
  S[3] = (quantile(obsData, 0.75) + quantile(obsData, 0.25) - 2*quantile(obsData, 0.5)) / S[2]
  
  S[4] = (quantile(obsData, 0.875) - quantile(obsData, 0.625) + quantile(obsData, 0.375) - quantile(obsData, 0.125)) / S[2]
  
  S[5] = S[1] * S[2]
  S[6] = S[1] * S[3]
  S[7] = S[1] * S[4]
  S[8] = S[2] * S[3]
  S[9] = S[2] * S[4]
  S[10] = S[3] * S[4]
  
  S[11:15] <- runif(5)
  
  return(S)
}

##########--------- Simulating data from g-and-k distribution for ABC ------------------################

# True Parameters
A = 3
B = 4
g = 2
k = 1

# Prior bounds
bounds <- cbind(c(0,0,0,0), c(10,10,10,10))

# Settings
nData = 10000 # Number of data points for each parameter value
nSamples = 2000 # Number of simulations from model
nParam = 4
nStats <- 15
nReal <- 100 # Number of times we repeat the experiment for each n_sim

# Initializing arrays
theta = array(0, dim = c(nReal, nSamples, nParam))
s_obs = array(0, dim = c(nReal, nStats))
sumStats = array(0, dim = c(nReal, nSamples, nStats))

for (i in 1:nReal) {
  set.seed(i)
  
  obsData_gk <- rgk(nData, A, B, g, k) # Sampling from g-and-k distribution
  s_obs[i,] <- gkStatistics(obsData_gk, nStats) # Computing observed statistics
  
  theta[i,,] <- array(runif(nSamples * nParam, 0, 10), dim = c(nSamples, nParam)) # Sampling parameter values from prior
  
  for (j in 1:nSamples) {
    data <- rgk(nData, theta[i,j,1], theta[i,j,2], theta[i,j,3], theta[i,j,4])
    sumStats[i,j,] <- gkStatistics(data, nStats)
  }
  print(i)
}

abc::abc(t(s_obs), theta, sumStats, tol, method = "loclinear", transf = 'logit', logit.bounds = bounds)$adj.values

saveRDS(s_obs, file = "gk_sobs_2000")

saveRDS(theta, file = "gk_param_2000")

saveRDS(sumStats, file = "gk_ssim_2000")




# generating reference data-set for gk

nSamples <- 10000
obsData_gk <- rgk(nData, A, B, g, k) # Sampling from g-and-k distribution
s_obs <- apply(readRDS("gk_sobs_2000"), 2, mean)

theta <- array(runif(nSamples * nParam, 0, 10), dim = c(nSamples, nParam))
sumStats = array(0, dim = c( nSamples, 15))
for (j in 1:nSamples) {
  data <- rgk(nData, theta[j,1], theta[j,2], theta[j,3], theta[j,4])
  sumStats[j,] <- gkStatistics(data, nStats)
  print(j)
}

true_samples <- abc::abc(s_obs[1:4], theta, sumStats[,1:4], tol=0.01, method = "loclinear", transf = 'logit', logit.bounds = bounds)$adj.values 

apply(true_samples, 2, mean)
saveRDS(true_samples, "gk_trueSamples")










